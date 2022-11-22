
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

filesD <- 'resultsWGCNA/'
inputData <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
inputNet <- paste0(filesD, 'network_manual_construction.RData')
load(inputData)
load(inputNet)

plotTOMF <- paste0(filesD, 'TOM heatmap.pdf')
sft_power <- 12
enhance_power <- 5

#Visualize gene network
dissTOM <- 1-TOMsimilarityFromExpr(datExpr, power = sft_power)
# Transform dissTOM with a power to make moderately strong connections more visible
plotTOM <- dissTOM^enhance_power
# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA

tiff(file = plotTOMF)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
invisible(dev.off())
save(plotTOM, file = paste0(filesD, 'temp.RData'))
#Vizualize eigengene network
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
condition <- as.data.frame(sampleTable$condition)
names(condition) <- 'condition'

MET <- orderMEs(MEs, condition)
par(cex = 0.9)
plotEigengeneNetworks(MET, '', marDendro= c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                      cex.lab= 0.8, xLabelsAngle = 90)