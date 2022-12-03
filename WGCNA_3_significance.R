#Weighted Gene Coexpression Network Analysis part 3, to relate modules to external
#clinical traits and gene significance
#Based on Peter Langfelder tutorials

suppressPackageStartupMessages({
  library(rstudioapi, quietly=T)
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
  #library(org.Hs.eg.db, quietly = T)
  library(anRichment)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

filesD <- 'results_WGCNA/'
inputData <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
inputNet <- paste0(filesD, 'network_manual_construction.RData')
load(inputData)
load(inputNet)

modtraitF <- paste0(filesD, "module-trait_relation.jpeg")
intramodFpref <- paste0(filesD, 'MMandGS/intramodular_')
geneInfoF <- paste0(filesD, 'genes_info.tsv')
finalresR <- paste0(filesD, 'summary_sigMods.RData')
gNumberF <- paste0(filesD, 'genes_in_significant_modules.txt')
allGenesListF <- paste0(filesD, 'genes_in_modules/all_modules_entrez.txt')

#Parameters
nSamples <- nrow(datExpr)
alpha <- 0.05

#A) Quantify module trait associations

#Transform qualitative traits to binary code
levels(sampleTable$condition) <- c(1,0)
colnames(sampleTable)[4] <- 'is.Affected'

#We recode the other categorical covariates to check if there is correlation
#It shouldnt because we already took them into account
levels(sampleTable$gender) <- c(1,0)
colnames(sampleTable)[1] <- 'is.Female'
levels(sampleTable$age) <- c(1,2,0)
levels(sampleTable)[2] <- 'age_stage'
#Probar añadiendo una columna con las edades reales

#Calculate correlations and p values taking into account # observations
moduleTrait <- corAndPvalue(MEs, sampleTable)
moduleTraitCor <- moduleTrait$cor
moduleTraitPvalue <- moduleTrait$p

#Get all together:
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #give shape

#Graphical representation:
jpeg(file = modtraitF, width = 700, height = 1600, quality = 100)
par(mar = c(6, 8.5, 3, 3))
title <- "Module-trait relationships"
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sampleTable),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, main = title)
invisible(dev.off())

#B) Focus on trait of interest: condition
condition <- as.data.frame(sampleTable$is.Affected)
names(condition) <- 'condition'
moduleNames <- substring(names(MEs), 3) #3 to take away prefix

#correlation of each gene in each sample to each module
modMemb <- corAndPvalue(datExpr, MEs)
geneModuleMembership <- modMemb$cor
geneModuleMembership <- as.data.frame(geneModuleMembership)
names(geneModuleMembership) <- paste("MM", moduleNames, sep="_")

MMPvalue <- modMemb$p
MMPvalue <- as.data.frame(MMPvalue)
names(MMPvalue) <- paste("p_MM", moduleNames, sep="_")

#correlation of each gene to each condition, thus, gene significance
geneSignif <- corAndPvalue(datExpr, condition)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(condition), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(condition), sep="")

#C) Intramodular analysis

#First get significant modules
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$is.Affected < alpha)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFpref, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())
}

#D) Summary output
#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste(filesD,"genes_in_modules/genes-", mod, ".txt", sep="");
  write.table(as.data.frame(modGenes), file = filename, quote = F,
              row.names = FALSE, col.names = FALSE)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)
}
#Gene information data frame
genes <- names(datExpr)
symbol <- mapIds(org.Hs.eg.db, keys = genes, column = 'SYMBOL', 
                 keytype = 'ENSEMBL', multiVals = 'first')
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                        geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, condition, method = "p")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

# Order the genes in the geneInfo variable first by module color, 
#then by geneTraitSignificance
geneOrder <- order(geneInfo$moduleColor, -abs(geneInfo$GS.condition))
geneInfo <- geneInfo[geneOrder, ]

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)


###ANNOTATION
#Save the genes present in each significant module, but this time in EntrezID
entrez <- mapIds(org.Hs.eg.db, keys = genes, column = 'ENTREZID', keytype = 'ENSEMBL', 
                 multiVals = 'first')
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- entrez[modGenes]
  filename <- paste(filesD,"genes_in_modules/genes-Entrez-", mod, ".txt", sep="");
  write.table(as.data.frame(modGenes), file = filename, quote = F,
              row.names = FALSE, col.names = FALSE)
}
#Save the rest of the genes as background for the analysis
write.table(as.data.frame(entrez), file = allGenesListF,quote = F, row.names=F,
            col.names=F)

#Enrichment GO analysis using WGCNA functions
GOcollection <- buildGOcollection(organism = "human")
for (mod in sigMods){
  active <- entrez[moduleColors==mod]
  enrichmnt <- enrichmentAnalysis(active = active, inactive = entrez,
                                           refCollection = GOcollection,
                                           useBackground = "intersection",
                                           threshold = 1e-4,
                                           thresholdType = "Bonferroni");
  enr_table <- enrichmnt$enrichmentTable
  terms <- enr_table$dataSetName
  filename <- paste0(filesD,"enriched_terms/terms_", mod, ".txt")
  write.table(as.data.frame(terms), file = filename, quote = F,
              row.names = FALSE, col.names = FALSE)
}

#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)
