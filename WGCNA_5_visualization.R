
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

plotTOMF <- paste0(filesD, 'TOM heatmap.tiff')
dendroEigenF <- paste0(filesD, 'Eigengene network dendrogram.tiff')
heatEigenF <- paste0(filesD, 'Eigengene network heatmap.tiff')

sft_power <- 12
enhance_power <- 5

#Visualize gene network
dissTOM <- 1-TOMsimilarityFromExpr(datExpr, power = sft_power)
# Transform dissTOM with a power to make moderately strong connections more visible
plotTOM <- dissTOM^enhance_power
# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA

tiff(file = plotTOMF)
title <- "Network heatmap plot, all genes"
TOMplot(plotTOM, geneTree, moduleColors, main = title)
invisible(dev.off())
#save(plotTOM, file = paste0(filesD, 'temp.RData'))

#Vizualize eigengene network
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
condition <- as.data.frame(sampleTable$condition)
names(condition) <- 'condition'
df <- cbind(MEs, condition)
levels(df$condition) <- c(1,0) #categorical to binary

MET <- orderMEs(df)
#Eigengene network dendrogram
tiff(filename = dendroEigenF)
title <- "Eigengene dendrogram"
plotEigengeneNetworks(MET, title, marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
invisible(dev.off())

#Eigengene network heatmap
tiff(filename = heatEigenF)
title <- "Eigengene adjacency heatmap"
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
invisible(dev.off())


