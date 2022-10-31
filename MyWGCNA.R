setwd("C:/Users/CBM/Desktop/RNAseq_bipolarDisorder")
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


#1. Load data
#Load variables
sampleTable <- read.table('configfile_pedlabels.txt', header=TRUE, row.names = 1, 
                          colClasses = rep('factor', 5))
sampleTable <- sampleTable[,c(3,4,5,2)]
names(sampleTable) #check condition, gender, age and PED are present

#Load expression data
df <- read.table(file = 'MyResults_DEG/counts_raw.tsv', sep = '\t', 
                 header = TRUE, row.names = 1)
#Filter low count genes, with less than 10 counts in 90% of features
keep <- rowSums(df > 10) > (ncol(df) * 0.9)
df <- df[keep,] #drastic drop from 28525 to 14938
#order as in metafile
df <- df %>%
  dplyr::select(rownames(sampleTable))
rm(keep)
#Get a DeSeq object
dds <- DESeqDataSetFromMatrix(countData = df, colData = sampleTable, 
                              design = ~ age + gender + PED + condition)

dds_norm <- vst(dds, blind = FALSE)
df_norm <- assay(dds_norm)
mm <- model.matrix(~condition, colData(dds_norm))
df_norm <- limma::removeBatchEffect(df_norm,
                                batch=dds_norm$PED, batch2=dds_norm$gender,
                                batch3=dds_norm$age, design=mm)

#Transpose for the WGCNA
datExpr <- as.data.frame(t(df_norm))

#Quality control
#check genes and samples with many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
#####
##cluster to detect outliers
sampleTree <- hclust(dist(datExpr), method = 'average')
sizeGrWindow(12,9)
pdf(file = "MyResults_WGCNA/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
invisible(dev.off())
#It doesn't seem to be any outlier
#####
#save relevant objects
#save(datExpr, sampleTable, file = 'MyResults_WGCNA/initialData_datExpr_sampleTable.RData')

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
# this line corresponds to using an R^2 cut-off of h
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
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft_s$fitIndices[,1], -sign(sft_s$fitIndices[,3])*sft_s$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft_s$fitIndices[,1], sft_s$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_s$fitIndices[,1], sft_s$fitIndices[,5], labels=powers, cex=cex1,col="red")
invisible(dev.off())

#Therefore we set the soft threshold parameter at 12 as it is the lowest at which
#we obtain a R^2 of 0.9
sft_power <- 12
minMod <- 15

#Build network
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
'moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "MyResults_WGCNA/netwrok_auto.RData")
'
###STEP BY STEP NETWORK CoNSTRUCTION
#coexpression similarity and adjacency
adjac <- adjacency(datExpr, power = sft_power)
#topological overlap matrix
tom <- TOMsimilarity(adjac, TOMType = 'signed', verbose = 5)
diss_tom <- 1-tom

geneTree <-hclust(as.dist(diss_tom), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = F, hang = 0.04)

#Cut the dendrogram
minMod <- 15
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

plot(METree, main = 'Clustering of Module Eigengenes')

#we set a height cut of 0.25 which corresponds to a correlation of 75%
MEDissThr <- 0.25
abline(h=MEDissThr, col = 'red')
merge <- mergeCloseModules(datExpr, dynColors, cutHeight = MEDissThr, verbose = 3)
mergedCols <- merge$colors
mergedMEs <- merge$newMEs


# Call the network topology analysis function


#3. Relating modules to external clinical traits and identifying important genes

#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
