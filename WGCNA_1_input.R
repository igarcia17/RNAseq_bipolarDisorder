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
