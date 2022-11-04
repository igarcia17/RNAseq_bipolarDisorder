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
# Call the network topology analysis function


#3. Relating modules to external clinical traits and identifying important genes

#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
