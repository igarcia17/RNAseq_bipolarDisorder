#Weightd Gene Coexpression Network Analysis part 2, for network construction
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

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

#Parameters
sign <- T #signed or unsigned network
Rcutoff <- 0.9 #R^2 cut off
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
minMod <- 15 #Set minimum module size

#Input/output
resD <- 'resultsWGCNA'
input <- paste0(resD,'initialData_datExpr_sampleTable.RData')
load(file = input)
unsignedSFTF <- paste0(resD,"softT_unsigned.tiff")
signedSFTF <- paste0(resD,"softT_signed.tiff")
autodendroF <- paste0(resD,"cluster_dendro_auto_min", minMod, ".tiff")

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

#PARAMETERS
#set the soft threshold parameter at 12 as it is the lowest at which
#we obtain a R^2 of 0.9
sft_power <- 12


#Build network in automatic way
#####
net <- blockwiseModules(datExpr, networkType = 'signed', minModuleSize = minMod,
                          TOMType = 'signed', power = 12, randomSeed = 1,
                          numericLabels = T, verbose = 5)

moduleColors <- labels2colors(net$colors)
geneTree <- net$dendrograms[[1]]
#Plot
tiff(file = autodendroF)
blocks <- moduleColors[net$blockGenes[[1]]]
plotDendroAndColors(geneTree, blocks,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
invisible(dev.off())
#save results
moduleLabels <- net$colors
MEs <- net$MEs
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "MyResults_WGCNA/network_auto.RData")
#####

###STEP BY STEP NETWORK CoNSTRUCTION
#coexpression similarity and adjacency
adjac <- adjacency(datExpr, power = sft_power)
#topological overlap matrix
tom <- TOMsimilarity(adjac, TOMType = 'signed', verbose = 5)
diss_tom <- 1-tom
geneTree <- hclust(as.dist(diss_tom), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = F, hang = 0.04)

#Save this objects, as running previous commands is very computationally expensive
save(tom, diss_tom, geneTree, file = 'MyResults_WGCNA/tom_dissTom_geneTree.RData')
load(file='MyResults_WGCNA/tom_dissTom_geneTree.RData')

#Cut the dendrogram
dynamicMods <- cutreeDynamic(geneTree, distM = diss_tom, deepSplit =2, 
                             pamRespectsDendro =F, minClusterSize = minMod)
table(dynamicMods)

dynColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynColors, 'Dynamic Tree Cut', dendroLabels = F, 
                    hang=0.03, addGuide = T, guideHang = 0.05, 
                    main = 'Gene dendrogram and module colors')

#Merge modules with similar expression profile
MEList <- moduleEigengenes(datExpr, colors = dynColors)
MEs <-MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = 'average')
#we set a height cut of 0.25 which corresponds to a correlation of 75%
MEDissThr <- 0.25

pdf(file = 'MyResults_WGCNA/Clustering_of_modules.pdf', width = 20, height = 8)
plot(METree, main = 'Clustering of Module Eigengenes')
abline(h=MEDissThr, col = 'red')
invisible(dev.off())

#Merge modules with high correlation among them
merge <- mergeCloseModules(datExpr, dynColors, cutHeight = MEDissThr, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
#Plot the modules with and without the merging
pdf(file='MyResults_WGCNA/Gene_dendrogram.pdf', wi = 20, he = 6)
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
save(MEs, moduleLabels, moduleColors, geneTree, file=
       'MyResults_WGCNA/network_manual_construction.RData')
