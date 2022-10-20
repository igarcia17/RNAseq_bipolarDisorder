setwd("C:/Users/CBM/Desktop/RNAseq_bipolarDisorder")

library(WGCNA)
options(stringsAsFactors = FALSE)

#1. Load data
df <- read.table(file = 'MyResults_DEG/counts_normalized.tsv', sep = '\t', header = FALSE)


#2. Network construction and module detection  

#3. Relating modules to external clinical traits and identifying important genes

#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
