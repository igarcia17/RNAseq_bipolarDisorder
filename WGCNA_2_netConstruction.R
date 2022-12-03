#Weighted Gene Coexpression Network Analysis part 2, for network construction
#Multiple methods are shown, though at the end the network of interest is manually
#built, soft threshold and signed. Bsed on the tutorials by Peter Langfelder

#Part 2: Network construction and module detection
suppressPackageStartupMessages({
  library(rstudioapi, quietly=T)
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
})
options(stringsAsFactors = F)
enableWGCNAThreads()

#Parameters
sign <- T #signed or unsigned network
Rcutoff <- 0.9 #R^2 cut off
powers <- c(c(1:10), seq(from = 12, to=20, by=2)) #possible thresholds
minMod <- 15 #Set minimum module size for preliminar-manual and automatic construction
autoNC <- F #make automatic network construction
MEDissThr <- 0.25 #threshold to which merge modules in manual construction
#height cut of 0.25 which corresponds to a correlation of 75% among modules

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

#Input/output
filesD <- 'resultsWGCNA/'
input <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
load(file = input)

unsignedSFTF <- paste0(filesD,"softT_unsigned.tiff")
signedSFTF <- paste0(filesD,"softT_signed.tiff")
autodendroF <- paste0(filesD,"cluster_dendro_auto_min", minMod, ".tiff")
autoR <- paste0(filesD, "network_auto.RData")
tempRdata <- paste0(filesD, 'tom_dissTom_geneTree.RData')
geneTree1F <- paste0(filesD, 'Preliminar_gene_clustering.pdf')
clustermodsF <- paste0(filesD, 'Clustering_of_modules.pdf')
dendroF <- paste0(filesD, 'Gene_dendrogram.pdf')
resultsR <- paste0(filesD, 'network_manual_construction.RData')

#Parameters: inspection of scale free topology by network type
if (sign){
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed', 
                           verbose = 5)
  filename <- signedSFTF
} else {
  sft <- pickSoftThreshold(datExpr, powerVector = powers, 
                           networkType = 'unsigned', verbose = 5)
  filename <- unsignedSFTF
}

#Plot how does data fit the scale free network model and mean connectivity
#based on the threshold used.
x <- sft$fitIndices[,1]
xlab <- "Soft Threshold (power)"
y1 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
ylab1 <- "Scale Free Topology Model Fit, signed R^2"
y2 <- sft$fitIndices[,5]
ylab2 <- "Mean Connectivity"

cex1 <- 0.9
tiff(file = filename, units="in", width=5, height=5, res=300)
par(mfrow = c(1,2))

plot(x, y1, xlab= xlab, ylab= ylab1, type="n", main = paste("Scale independence"))
text(x, y1, labels=powers, cex=cex1, col="red")
abline(h=Rcutoff, col="red")

plot(x, y2, xlab=xlab, ylab=ylab2, type="n", main = ylab2)
text(x, y2, labels=powers, cex=cex1, col="red")
invisible(dev.off())

#Set the soft threshold parameter at 12 as it is the lowest at which
#we obtain a R^2 of 0.9
sft_power <- 12

#Build network in automatic way
if (autoNC){
  net <- blockwiseModules(datExpr, networkType = 'signed', minModuleSize = minMod,
                          power = sft_power, randomSeed = 1,
                          numericLabels = F, verbose = 5)
  moduleColors <- labels2colors(net$colors)
  geneTree <- net$dendrograms[[1]]
  
#Plot
  tiff(file = autodendroF)
  blocks <- moduleColors[net$blockGenes[[1]]]
  plotDendroAndColors(geneTree, blocks, "Module colors",dendroLabels = FALSE, 
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  invisible(dev.off())
  
#save results
  moduleLabels <- net$colors
  MEs <- net$MEs
  save(MEs, moduleLabels, moduleColors, geneTree, file = autoR)
}

##STEP BY STEP NETWORK CONSTRUCTION

#co expression similarity and adjacency
adjac <- adjacency(datExpr, power = sft_power)
#topological overlap matrix: compute how much directly and indirectly two genes are
tom <- TOMsimilarity(adjac, TOMType = 'signed', verbose = 5)
diss_tom <- 1-tom
#Gene tree plot
geneTree <- hclust(as.dist(diss_tom), method = 'average')
title <- 'Gene clustering on TOM-based dissimilarity'
plot(geneTree, xlab = '', sub = '', main = title, labels = F, hang = 0.04)

#Save this objects, as running previous commands is very computationally expensive
#save(tom, diss_tom, geneTree, file = tempRdata)
#load(file=tempRdata)

#Cut the dendrogram
dynamicMods <- cutreeDynamic(geneTree, distM = diss_tom, deepSplit =2, 
                             pamRespectsDendro =F, minClusterSize = minMod)
table(dynamicMods)
dynColors <- labels2colors(dynamicMods)
title <- 'Gene dendrogram and module colors'
pdf(geneTree1F)
plotDendroAndColors(geneTree, dynColors, 'Dynamic Tree Cut', dendroLabels = F, 
                    hang=0.03, addGuide = T, guideHang = 0.05, main = title)
invisible(dev.off())
#Merge modules with similar expression profile
MEList <- moduleEigengenes(datExpr, colors = dynColors)
MEs <-MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = 'average')

#Plot
title <- 'Clustering of Module Eigengenes'
pdf(file = clustermodsF, width = 20, height = 8)
plot(METree, main = title)
abline(h=MEDissThr, col = 'red')
invisible(dev.off())

#Merge modules with high correlation among them
merge <- mergeCloseModules(datExpr, dynColors, cutHeight = MEDissThr, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

#Plot the modules with and without the merging
pdf(file=dendroF, wi = 20, he = 6)
plotDendroAndColors(geneTree, cbind(dynColors, mergedColors), 
                    c('Dynamic Tree Cut','Merged dynamic'),
                    dendroLabels = F, hang = 0.03, guideHang = 0.05, addGuide = T)
invisible(dev.off())

#Save the results
moduleColors <- mergedColors
colorOrder <- c('grey', standardColors(n=NULL))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

#Save module colors and labels for subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file=resultsR)
