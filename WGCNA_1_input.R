#Weighted Gene Coexpression Network Analysis first part, for data loading
#Based on tutorials of WGCNA package author Peter Langfelder

#Part 1: Data loading
suppressPackageStartupMessages({
  library(rstudioapi)
  library(dplyr, quietly = T)
  library(DESeq2, quietly = T)
  library(WGCNA, quietly = T)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
covs <- T

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

#Input
inputC <- 'configfile.txt'
inputTSV <- 'results_DGE/counts_raw.tsv'

#Outputs
resD <- 'results_WGCNA/'
samtreeF <- paste0(resD,"sampleClustering.jpeg")
finalresults <- paste0(resD,'initialData_datExpr_sampleTable.RData')

#Load covariates
sampleTable <- read.table(inputC, header=TRUE, row.names = 1, colClasses = 
                            rep('factor', 5))
sampleTable <- sampleTable[,c(3,4,5,2)]
#names(sampleTable)

#Load expression data
df <- read.table(file = inputTSV, sep = '\t', header = TRUE, row.names = 1)

#Filter low count genes, with less than 10 counts in 90% of features
keep <- rowSums(df > 10) > (ncol(df) * 0.9)
df <- df[keep,] #drastic drop from 28525 to 14938

#order as in configfile
df <- df %>%
  dplyr::select(rownames(sampleTable))
rm(keep) #Remove because big size

#Get a DeSeq object
dds <- DESeqDataSetFromMatrix(countData = df, colData = sampleTable, 
                              design = ~ age + gender + PED + condition)
#Remove effect of covariates if covs is TRUE
if(covs){
dds_norm <- vst(dds, blind = FALSE)
df_norm <- assay(dds_norm)
mm <- model.matrix(~condition, colData(dds_norm))
df_norm <- limma::removeBatchEffect(df_norm,
                                    batch=dds_norm$PED, batch2=dds_norm$gender,
                                    batch3=dds_norm$age, design=mm)
datExpr <- as.data.frame(t(df_norm))
} else {
  dds_norm <- vst(dds, blind = FALSE)
  df_norm <- assay(dds_norm)
  datExpr <- as.data.frame(t(df_norm))
}

#Quality control
#Remove genes and samples with many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK){
  # Print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    remG <- paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")
    printFlush(paste0("Removing genes: ", remG))
  if (sum(!gsg$goodSamples)>0)
    remS <- paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")
    printFlush(paste0("Removing samples: ", remS));
  # Remove the bad genes and samples from the data:
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

#Cluster to detect outliers in samples
sampleTree <- hclust(dist(datExpr), method = 'average')
jpeg(file = samtreeF, width = 700, height = 700, quality = 300)
par(mar = c(0,4,2,0))
title <- "Sample clustering to detect outliers"
plot(sampleTree, main = title, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
invisible(dev.off())
#It doesn't seem to be any outlier

#save relevant objects
save(datExpr, sampleTable, file = finalresults)

