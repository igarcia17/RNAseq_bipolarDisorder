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
  library(org.Hs.eg.db)
  library(fgsea, quietly = T)
  library(msigdbr, quietly = T)
})

data <- read.delim('MyResults_DEG/all_genes.tsv', sep= "\t", header=T, row.names = 1)

#Calculate GSEA

dat <- data$stat #define why use stat
names(dat) <- as.character(rownames(data))
dat <- sort(dat, decreasing=TRUE)

#Calculate GSEA and write tables of results

category <- 'H'
db_sets <- msigdbr(species = 'Homo sapiens', category = category)%>% 
  dplyr::select(gs_name, ensembl_gene)
head(db_sets) #each gene associated with each msig group

set.seed(1)
egs <- GSEA(geneList = dat, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets)
head(egs@result)
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

write.table(egs_df, file = 'MyResults_GSEA/GSEA_results.txt', sep= "\t", 
            quote = F, row.names = T)

##Dotplot / BarPlot

jpeg(file = "MyResults_GSEA/Dotplot.jpeg", units = 'in', width = 15, height = 10, 
     res = 300)
par(mar = c(2, 2, 2, 5)) 
dotplot(egs, x = "GeneRatio", color = "p.adjust", showCategory = 20, 
            font.size = 15)
invisible(dev.off())

##Gene-concept network

jpeg(file = "MyResults_GSEA/Gene_concept_net.jpeg", units = 'in', 
     width = 15, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
cnetplot(egs, categorySize="p.adjust", font.size = 15, colorEdge = T)
invisible(dev.off())

##Ridgeline plot

jpeg(file = "MyResults_GSEA/GSEA_ridge.jpeg", units = 'in', width = 15, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
ridgeplot(egs, fill="p.adjust", orderBy= 'NES', core_enrichment = T)
invisible(dev.off())

##Upset plot (of the 20 first terms)

genes_20first <- as.data.frame(as.factor(head(egs@result$core_enrichment, 20)))
lista_20first <- list()
for (i in 1:nrow(genes_20first)){
    lista_20first[[i]] <- unlist(strsplit(as.character(genes_20first[i,1]),split="/"))   
}
uniq_genes <- as.character(unique(names(dat))) #es necesario?

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


jpeg(file = "MyResults_GSEA/GSEA_upset_20first.jpeg", units = 'in', width = 15, height = 10, res = 300)
    upset(mat_20first, nsets=10, order.by="freq", sets.bar.color="skyblue")
invisible(dev.off())

##################################

###Plot the GSEA if terms are provided and if not, plot the first 5 more abundant terms

for (j in 1:5){
  pl <- gseaplot2(egs, geneSetID=j, title = egs$Description[j], base_size=40, color="red")
  desc <- gsub(" ", "_", egs$Description[j], fixed = TRUE) 
  filename <- paste('MyResults_GSEA/',desc, "_GSEA.jpeg", sep ="")
  ggsave(pl, file=filename, device = "jpeg", units= "in", height = 15, width = 20)
}

gseap <- gseaplot2(egs, geneSetID = 'todos', pvalue_table = T)
ggsave(gseap, file='All gseaplots', device = "jpeg", units= "in", height = 15, width = 20)
