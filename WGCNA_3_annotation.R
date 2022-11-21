#Weighted Gene Coexpression Network Analysis part 3, to relate modules to external
#clinical traits and driver gene identification and annotation
#Based on Peter Langfelder tutorials

suppressPackageStartupMessages({
  library(rstudioapi, quietly=T)
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
  library(org.Hs.eg.db, quietly = T)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

filesD <- 'resultsWGCNA/'
inputData <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
inputNet <- paste0(filesD, 'network_manual_construction.RData')
load(inputData)
load(inputNet)

modtraitF <- paste0(filesD, "Module-trait_relation.pdf")
intramodFpref <- paste0(filesD, 'intramodular_')

#Parameters
nSamples <- nrow(datExpr)
alpha <- 0.01

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
cp <- corAndPvalue(MEs, sampleTable)
moduleTraitCor <- cp$cor
moduleTraitPvalue <- cp$p

#Get all together:
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor) #give shape

#Graphical representation:
pdf(file = modtraitF, height = 20, width = 10)
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
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p")) #not the same as moduleTraitCor
names(geneModuleMembership) <- paste("MM", moduleNames, sep="_")

#p value of last correlation
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(MMPvalue) <- paste("p_MM", moduleNames, sep="_")

#correlation of each gene to each condition
geneTraitSignificance <- as.data.frame(cor(datExpr, condition, use = "p"))
names(geneTraitSignificance) <- paste("GS.", names(condition), sep="")

#p values of last correlation
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(GSPvalue) <- paste("p.GS.", names(condition), sep="")

#C) Intramodular analysis

#First get significant modules
modTraitP_df <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(modTraitP_df, modTraitP_df$is.Affected < alpha)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (i in sigMods){
  mod <- i
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

#3.4 Summary output
#Gene information data frame
genes <- names(datExpr)
symbol <- mapIds(org.Hs.eg.db, keys = genes, column = 'SYMBOL', 
                 keytype = 'ENSEMBL', multiVals = 'first')
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                        geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
modOrder <- order(-abs(cor(MEs, condition, use = "p")))

#Relate to DEGs


#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
