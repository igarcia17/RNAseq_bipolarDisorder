setwd("C:/Users/Asus/OneDrive/Escritorio/RNAseq_bipolarDisorder")
suppressPackageStartupMessages({
library(dplyr, quietly = T)
library(DESeq2, quietly = T)
library(magrittr, quietly = T)
library(ggplot2, quietly = T)
library(WGCNA, quietly = T)
#library(limma, quietly = T)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#2. Network construction and module detection
load(file = 'MyResults_WGCNA/initialData_datExpr_sampleTable.RData')

#Determine parameters: soft threshold
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
#####
#For an unsigned network:
sft_un <- pickSoftThreshold(datExpr, powerVector = powers, 
                            networkType = 'unsigned', verbose = 5)
tiff(file = "MyResults_WGCNA/softT_unsigned.tiff", units="in", width=5, 
     height=5, res=300)
par(mfrow = c(1,2))
cex1 <- 0.9
plot(sft_un$fitIndices[,1], -sign(sft_un$fitIndices[,3])*sft_un$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft_un$fitIndices[,1], -sign(sft_un$fitIndices[,3])*sft_un$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
plot(sft_un$fitIndices[,1], sft_un$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_un$fitIndices[,1], sft_un$fitIndices[,5], labels=powers, cex=cex1,col="red")
invisible(dev.off())
#####
#For a signed network:
sft_s <- pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed', verbose = 5)

tiff(file = "MyResults_WGCNA/softT_signed.tiff")
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft_s$fitIndices[,1], -sign(sft_s$fitIndices[,3])*sft_s$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))
text(sft_s$fitIndices[,1], -sign(sft_s$fitIndices[,3])*sft_s$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft_s$fitIndices[,1], sft_s$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_s$fitIndices[,1], sft_s$fitIndices[,5], labels=powers, cex=cex1,col="red")
invisible(dev.off())
#####
#PARAMETERS
#set the soft threshold parameter at 12 as it is the lowest at which
#we obtain a R^2 of 0.9
sft_power <- 12
#Set minimum module size
minMod <- 15

#Build network in automatic way
#####
net <- blockwiseModules(datExpr, networkType = 'signed', minModuleSize = minMod,
                          TOMType = 'signed', power = 12, randomSeed = 1,
                          numericLabels = T, verbose = 5)
#Plot
tiff(file = "MyResults_WGCNA/cluster_dendo_auto_min15.tiff")
mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
invisible(dev.off())
#save results
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
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
