setwd("C:/Users/Asus/OneDrive/Escritorio/RNAseq_bipolarDisorder")
suppressPackageStartupMessages({
  library(dplyr, quietly = T)
  library(magrittr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
load('Myresults_WGCNA/initialData_datExpr_sampleTable.RData')
load('MyResults_WGCNA/network_manual_construction.RData')
#3. Relating modules to external clinical traits and identifying important genes
#Quantify module trait associations
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

levels(sampleTable$condition) <- c(1,0)
colnames(sampleTable)[4] <- 'is.Affected'
#We recode the other categorical covariates to check if there is correlation
#It shouldnt because we already took them into account
levels(sampleTable$gender) <- c(1,0)
colnames(sampleTable)[1] <- 'is.Female'
levels(sampleTable$age) <- c(1,2,0)
levels(sampleTable)[2] <- 'age_stage'
#Probar añadiendo una columna con las edades reales
moduleTraitCor <- cor(MEs, sampleTable, use = "p")
#As previously thought, the cor between covariates is very small
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
#Como sacar el p valor ajustado

cp = corAndPvalue(MEs, sampleTable)
moduleTraitCor2 = cp$cor
moduleTraitPvalue2 = cp$p

# Call the network topology analysis function


#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
