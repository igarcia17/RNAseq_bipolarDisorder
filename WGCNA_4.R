#Weighted Gene Coexpression Network Analysis part 4, for identification
#of driver genes and comparison with previous DGE results.

suppressPackageStartupMessages({
  library(rstudioapi, quietly=T)
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
})
options(stringsAsFactors = F)
enableWGCNAThreads()

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

filesD <- 'results_WGCNA/'
inputData <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
inputNet <- paste0(filesD, 'network_manual_construction.RData')
inputSigMods <- paste0(filesD, 'summary_sigMods.RData')
inputSigGenes <- 'results_DGE/0.05_sig_padj.tsv'
load(inputData)
load(inputNet)
load(inputSigMods)

vennF <- paste0(filesD, 'common_modules_fromWGCNA_DEGs.tiff')

#Relate to previous DEG results
resDGE <- read.table(file = inputSigGenes, sep = '\t', header = TRUE, row.names=1)
resDGE <- resDGE[order(row.names(resDGE)),]
DEGs <- rownames(resDGE)

DEGinfo <- filter(geneInfo, rownames(geneInfo) %in% DEGs)
DEGinfo <- DEGinfo[order(row.names(DEGinfo)),]

#Why are there 4 genes missing?
missingGenes <- setdiff(DEGs, rownames(DEGinfo))
missGenesDGEs <- filter(resDGE, rownames(resDGE) %in% missingGenes)
#Looking at baseMean, we can see that they have probably been filtered out for the 
#WGCNA because of low counts, or because they are not significant enough

DEGmods_list <- unique(DEGinfo$moduleColor)
DEGmods_table <- table(DEGinfo$moduleColor)
#What percentage of DEGS are in a significant module?

#Pie chart, DEGs in module out of sigMods are in separate cathegory

#What significant modules appear in the DEGs? Venn diagram
library(VennDiagram)
venn.diagram(x = list(DEGmods_list, sigMods), category.names = 
               c('Modules of significant genes', 'Significant modules'),
             filename = vennF)
USING GGVEN PACKAGE
https://stackoverflow.com/questions/25019794/venn-diagram-with-item-labels