setwd("C:/Users/CBM/Desktop/RNAseq_bipolarDisorder")
library(dplyr, quietly = T)
library(DESeq2, quietly = T)
library(magrittr, quietly = T)
library(ggplot2, quietly = T)
library(WGCNA, quietly = T)
#library(limma, quietly = T)
options(stringsAsFactors = FALSE)

#1. Load data
#Load variables
sampleTable <- read.table('configfile_pedlabels.txt', header=TRUE, row.names = 1, 
                          colClasses = c('character','factor', 'factor', 
                                         'factor', 'factor'))
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
#Remove batch effect of covariates
mm <- model.matrix(~condition, data = sampleTable)
df <- limma::removeBatchEffect(df,
                                batch=sampleTable$PED, batch2=sampleTable$gender,
                                batch3=sampleTable$age, design=mm)

#Transpose for the WGCNA
datExpr <- as.data.frame(t(df))

#Quality control
#check genes and samples with many missing values
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

#cluster to detect outliers
sampleTree <- hclust(dist(datExpr0), method = 'average')
sizeGrWindow(12,9)
#pdf(file = "MyResults_WGCNA/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#It doesn't seem to be any outlier






#save relevant objects
save(datExpr, sampleTable, file = 'MyResults_WGCNA/initialData_datExpr_sampleTable.RData')

#2. Network construction and module detection



#3. Relating modules to external clinical traits and identifying important genes

#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
