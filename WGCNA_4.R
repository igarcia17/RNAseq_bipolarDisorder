
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
inputSig <- paste0(filesD, 'summary_sigMods.RData')
load(inputData)
load(inputNet)
load(inputSig)


#Relate to DEGs

#Load significant genes, make it a list
#relate the significant genes to their associated color
#Graph, like a pie



