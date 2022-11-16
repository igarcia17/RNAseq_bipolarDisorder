setwd("C:/Users/Asus/OneDrive/Escritorio/RNAseq_bipolarDisorder")
suppressPackageStartupMessages({
  library(dplyr, quietly = T)
  library(magrittr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
load('MyResults_WGCNA/initialData_datExpr_sampleTable.RData')
load('MyResults_WGCNA/network_manual_construction.RData')

nSamples <- nrow(datExpr)
alpha <- 0.01

#3. Relating modules to external clinical traits and identifying important genes
#3.1 Quantify module trait associations
#Transform qualitative traits to binary code
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
#It gets the same results as the following code, which takes
#into account the number of tests, thus, the adjusted p value:
#cp <- corAndPvalue(MEs, sampleTable)
#moduleTraitCor2 <- cp$cor
#moduleTraitPvalue2 <- cp$p

#Get all together:
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor) #give shape

#####
#Graphical representation:

pdf(file = "MyResults_WGCNA/Module-trait_relation.pdf",
     height = 20, width = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sampleTable),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               main = paste("Module-trait relationships"))
invisible(dev.off())
#####
#3.2 Focus on trait of interest: condition
condition <- as.data.frame(sampleTable$is.Affected)
names(condition) <- 'condition'
moduleNames <- substring(names(MEs), 3) #3 to take away prefix

#correlation of each gene in each sample to each module
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p")) #not the same as moduleTraitCor
names(geneModuleMembership) <- paste("MM", moduleNames, sep="_")
#p value of last correlation
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(MMPvalue) <- paste("p_MM", moduleNames, sep="_")

#correlation of each gene to each condition
geneTraitSignificance <- as.data.frame(cor(datExpr, condition, use = "p"))
names(geneTraitSignificance) <- paste("GS.", names(condition), sep="")
#p values of last correlation
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(GSPvalue) <- paste("p.GS.", names(condition), sep="")

#3.3 Intramodular analysis
#First get significant modules
modTraitP_df <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(modTraitP_df, modTraitP_df$is.Affected<alpha)
sigMods <- substring(rownames(sigMods), 3)
#Plot significance of each gene in every significant module
# against their correlation to the module.
for (i in sigMods){
  mod <- i
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  tiff(filename = paste('MyResults_WGCNA/intramodular_',mod,'.tiff', sep=''), units="in", 
       width=5, height=5, res=300)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", mod, "module"),
                     ylab = "Gene significance for condition",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())
}

#3.4 Summary output


#Relate to DEGs

# Call the network topology analysis function


#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
