#!/usr/bin/env Rscript
##
##R Script to make GSEA analysis
##The input must be a tab-separated file with 2 columns: first column must contain the gene IDs
##and the second one the log2FC values.
##
##THE ORDER OF THE COLUMNS IS VERY IMPORTANT!!
## 
## Eva Sacristán and Sandra González (GENGS CBMSO)

#####################
setwd("C:/Users/CBM/Desktop/RNAseq_bipolarDisorder")

suppressPackageStartupMessages({
  library(BiocManager, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(UpSetR, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(optparse, quietly = TRUE)
  library(org.Hs.eg.db)
  library(fgsea, quietly = T)
})

data <- read.table('MyResults_DEG/gsea_list.tsv', sep= "\t", header=F)
prefix <- 'GSEA'

##GENERATE PLOT
#####################

cat("\n Calculation GSEA and making some nice plots \n")
    
##Calculate GSEA

dat <- data$V2
names(dat) <- as.character(data$V1)
dat <- sort(dat, decreasing=TRUE)

#Calculate GSEA and write tables of results
set.seed(1)
egs <- gseGO(geneList = dat, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", 
             ont= "ALL")
egs_genename <- setReadable(egs, OrgDb = "org.Hs.eg.db")

write.table(egs, file = paste("tableGO_",prefix,".txt", sep =""), sep= "\t", quote = F)
write.table(egs_genename, file = paste("tableGO_",prefix,"_genename.txt", sep =""), sep= "\t", quote = F)

##Dotplot / BarPlot

jpeg(file = paste(prefix, "_dotplot.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    dotplot(egs, x = "GeneRatio", color = "pvalue", showCategory = 20, font.size = 15)
invisible(dev.off())

##Gene-concept network

jpeg(file = paste(prefix, "_gene_concept_net.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    cnetplot(egs_genename, categorySize="pvalue", foldChange=dat_sort, font.size = 15, colorEdge = T)
invisible(dev.off())

##Ridgeline plot

jpeg(file = paste(prefix, "_ridge.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    ridgeplot(egs, fill="pvalue")
invisible(dev.off())

##Heatplot

jpeg(file = paste(prefix, "_ridge.jpeg", sep =""), units = 'in', width = 20, height = 10, res = 300)
    heatplot(egs, foldChange=dat_sort)
invisible(dev.off())

##Upset plot (of the 20 first terms)

genes_20first <- as.data.frame(as.factor(head(egs@result$core_enrichment, 20)))
lista_20first <- list()
for (i in 1:nrow(genes_20first)){
    lista_20first[[i]] <- unlist(strsplit(as.character(genes_20first[i,1]),split="/"))   
}
uniq_genes <- as.character(unique(names(dat_sort)))

func_20first <- egs$Description[1:20]
mat <- matrix(0L, nrow = length(uniq_genes), ncol = length(func_20first)) 

for (i in 1:length(uniq_genes)) {
for (j in 1:length(func_20first)) {
gen <- uniq_genes[i]
if (gen %in% lista_20first[[j]]) {
mat[i,j] =  1
}}} 
 
mat_20first <- as.data.frame(mat)
colnames(mat_20first) <- func_20first
row.names(mat_20first) <- uniq_genes


jpeg(file = paste(prefix, "_upset_20first.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    upset(mat_20first, nsets=10, order.by="freq", sets.bar.color="skyblue")
invisible(dev.off())

##################################

###Plot the GSEA if terms are provided and if not, plot the first 5 more abundant terms

for (j in 1:5){
  pl <- gseaplot2(egs, geneSetID=j, title = egs$Description[j], base_size=40, color="red")
  desc <- gsub(" ", "_", egs$Description[j], fixed = TRUE) 
  filename <- paste(desc, "_", prefix, ".jpeg", sep ="")
  ggsave(pl, file=filename, device = "jpeg", units= "in", height = 15, width = 20)
}








## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = dat, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)
