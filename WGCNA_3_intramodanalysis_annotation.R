#Weighted Gene Coexpression Network Analysis part 3, to relate modules to external
#clinical traits and gene significance
#Based on Peter Langfelder tutorials
#Annotation based on:
#https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html#4_Over-Representation_Analysis_with_clusterProfiler_-_RNA-seq

suppressPackageStartupMessages({
  library(rstudioapi, quietly=T)
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
  library(org.Hs.eg.db, quietly = T)
  library(clusterProfiler)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#Path and files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

filesD <- 'results_WGCNA/'
inputData <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
inputNet <- paste0(filesD, 'network_manual_construction.RData')
load(inputData)
load(inputNet)

#For output
subD <- NULL
modtraitF <- paste0(filesD, subD, "module-trait_relation.jpeg")
intramodFprefix <- paste0(filesD, subD, 'MM_and_GS/intramodular_')
geneInfoF <- paste0(filesD, subD, 'genes_info.tsv')
finalresR <- paste0(filesD, subD, 'summary_sigMods.RData')
gNumberF <- paste0(filesD, subD, 'genes_in_significant_modules.txt')
allGenesListF <- paste0(filesD, subD, 'genes_in_modules/all_modules_entrez.txt')
annotationFprefix <- paste0(filesD, subD, 'annotation/enrich_results_')
genelistsFprefix <- paste0(filesD, subD, "genes_in_modules/genes-")

#Parameters
nSamples <- nrow(datExpr)
alpha <- 0.05
category <- 'C5'
subcategory <- NULL

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

#P adjusted?
#WGCNA alleviates the multiple testing problem inherent in microarray data analysis.
# Instead of relating thousands of genes to a microarray sample trait, it focuses on the relationship between a few (typically less than 10)
# modules and the sample trait. Toward this end, it calculates the eigengene significance (correlation between sample trait and eigengene)
# and the corresponding p-value for each module


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
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())
}

#D) Summary output
genes <- names(datExpr)
symbol <- mapIds(org.Hs.eg.db, keys = genes, column = 'SYMBOL', 
                 keytype = 'ENSEMBL', multiVals = 'first')

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(as.data.frame(modGenes), file = filename, quote = F,
              row.names = FALSE, col.names = FALSE)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)
}

#Gene information data frame
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
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

###ANNOTATION: Over Representation Analysis
#Make a data frame with each gene related to their assigned module
gene_module_df <- data.frame(rownames(geneInfo), geneInfo$moduleColor)
colnames(gene_module_df) <- c('gene', 'module')
#Retrieve data base
hs_msigdb_df <- msigdbr(species = 'Homo sapiens', category = category, 
                   subcategory = subcategory)%>% 
  dplyr::select(gs_name, ensembl_gene)
#Make a background gene set with all possible genes, to compare later on
background_set <- unique(as.character(gene_module_df$gene))

#For each significant module make a over representation analysis
set.seed(13)
for (mod in sigMods){
  #Make a list of genes in the module
  active_genes <- gene_module_df %>%
    dplyr::filter(module == mod) %>%
    dplyr::pull("gene")
  #Make analysis
  enriching <- enricher(gene = active_genes,
                        pvalueCutoff = 0.1,
                        pAdjustMethod = "BH",
                        universe = background_set,
                        TERM2GENE = hs_msigdb_df)
  enrich_res <- data.frame(enriching@result)
  #Save results
  write.table(enrich_res, file = paste0(annotationFprefix, mod, '.txt'),
              quote = F, row.names=F, col.names=T)
  enrich_res_sig <- enrich_res %>%
    dplyr::filter(p.adjust < alpha)
  #Just if there is any significant enriched term
  if(nrow(enrich_res_sig) != 0){
    write.table(enrich_res_sig, file = paste0(annotationFprefix, 'sig_', mod, '.txt'),
                quote = F, row.names=F, col.names=T)
    #Only make plots if there are more than 1 enriched terms
    if (nrow(enrich_res_sig) > 1){
      enrich_plot <- enrichplot::dotplot(enriching)
      upset_plot <- enrichplot::upsetplot(enriching)
      
      filename <- paste0(annotationFprefix, 'dotplot_', mod, '.jpeg')
      ggsave(enrich_plot, file=filename, device = "jpeg", units= "in", 
             height = 10, width = 15)
      
      filename <- paste0(annotationFprefix, 'upsetplot_', mod, '.jpeg')
      ggsave(upset_plot, file=filename, device = "jpeg", units= "in", 
             height = 15, width = 20)
    }
    } else {
      message <- paste0('The module ', mod, 
                        ' doesnt have significantly enriched terms')
      print(message)
    }
  #To check that process is going OK
  message <- paste0('Finished with module ', mod)
  print(message)
  }

