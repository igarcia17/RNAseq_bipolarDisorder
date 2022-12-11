## Based on the scripts by Eva Sacristán and Sandra
#González (GENGS CBMSO)

suppressPackageStartupMessages({
  library(rstudioapi)
  library(DESeq2)
  library(clusterProfiler, quietly = TRUE)
  library(msigdbr, quietly = T)
  library(UpSetR, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(ggplot2)
})

#Parameters
#Which database inside msigdbr?
category <- 'H'
subcategory <- NULL
#msigdbr_collections()
#Will you use stat parameter for ordering or the shrinked log2fold change?
statP <- F
#Plot the x top categories
topCat <- 5

#Paths
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
#Input
input <- 'results_DGE/deseq_objects.RData'

#Outputs
resD0 <- 'results_GSEA/'
if (statP){
  resD1 <- paste0(resD0,'stat/')
} else {
  resD1 <- paste0(resD0, 'log2fold/')
}
resD <- gsub(':','_',paste0(resD1,category,'_', subcategory, '/'))
if (!file.exists(resD)){
  dir.create(file.path(resD))
}

resTSV <- paste0(resD,'GSEA_results.txt')
dotplotF <- paste0(resD, "dotplot.jpeg")
geneconceptF <- paste0(resD,"gene_concept_net.jpeg")
ridgeF <- paste0(resD,"GSEA_ridge.jpeg")
upsetF <- paste0(resD,"upset_plot.jpeg")
gseaplotsF <- paste0(resD,'all_gseaplots.jpeg')

#1) Load data
load(input)

if(statP){
  res <- results(dds)
  res <- res[complete.cases(res),]
  dat <- res$stat
  names(dat) <- as.character(rownames(res))
  dat <- sort(dat, decreasing=TRUE)
} else {
  shrink <- lfcShrink(dds, coef = 12, type="apeglm", quiet =T)
  dat <- shrink$log2FoldChange
  names(dat) <- as.character(rownames(shrink))
  dat <- sort(dat, decreasing=TRUE)
}

#2) Calculate GSEA and write tables of results
#Get genes and categories
db_sets <- msigdbr(species = 'Homo sapiens', category = category, 
                   subcategory = subcategory)%>% 
  dplyr::select(gs_name, ensembl_gene)
head(db_sets) #each gene associated with each msig group

#Perform GSEA
set.seed(1)
egs <- GSEA(geneList = dat, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
#https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html 
#head(egs@result)
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

write.table(egs_df, file = resTSV, sep= "\t", quote = F, row.names = T)

#Reconsider the top category number if there are less terms than especified
if (dim(egs_df)[1] < topCat){
  topCat <- dim(egs_df)[1]
}

#Plot the results
##Dotplot
jpeg(file = dotplotF, units = 'in', width = 15, height = 10, 
     res = 300)
par(mar = c(2, 2, 2, 5)) 
title <- 'Dot plot with GSEA categories'
dotplot(egs, x = "GeneRatio", color = "p.adjust", showCategory = 20, 
            font.size = 15, title = title)
invisible(dev.off())

##Gene-concept network
jpeg(file = geneconceptF, units = 'in', width = 15, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
cnetplot(egs, categorySize="p.adjust", font.size = 15, colorEdge = T)
invisible(dev.off())

##Ridge line plot
jpeg(file = ridgeF, units = 'in', width = 15, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
ridgeplot(egs, fill="p.adjust", orderBy= 'NES', core_enrichment = T, 
          showCategory = topCat)
invisible(dev.off())

##Upset plot (of the 10 first terms)
#Save the genes in each category in a list
genes_top <- as.data.frame(as.factor(head(egs@result$core_enrichment, topCat)))
list_top <- list()
for (i in 1:topCat) {
    list_top[[i]] <- unlist(strsplit(as.character(genes_top[i,1]),split="/"))   
}

#Store all unique gene IDs
uniq <- as.character(unique(names(dat)))
#Get top functions
func_top <- egs$Description[1:topCat]
#Make sparse matrix with 1 for every gene in each category
mat <- matrix(0L, nrow = length(uniq), ncol = length(func_top)) 
for (gene in 1:length(uniq)) {
  for (func in 1:length(func_top)) {
    gen <- uniq[gene]
    if (gen %in% list_top[[func]]) {
      mat[gene,func] =  1
    }}} 

#Make a data frame
mat_df <- as.data.frame(mat)
colnames(mat_df) <- func_top
row.names(mat_df) <- uniq

#Plot
jpeg(file = upsetF, units = 'in', width = 15, height = 10, res = 300)
upset(mat_df, nsets=10, order.by="freq", sets.bar.color="skyblue")
invisible(dev.off())

###Plot the first 5 more abundant terms or all if there are less hits

for (j in 1:topCat){
  pl <- gseaplot2(egs, geneSetID=j, title = egs$Description[j], base_size=40, color="red")
  desc <- gsub(" ", "_", egs$Description[j], fixed = TRUE) 
  filename <- paste0(resD, desc, ".jpeg")
  ggsave(pl, file=filename, device = "jpeg", units= "in", height = 15, width = 20)
}

#For all top categories at once:
gseap <- gseaplot2(egs, geneSetID = 1:topCat, pvalue_table = F)
ggsave(gseap, file=gseaplotsF, device = "jpeg", units= "in", 
       height = 15, width = 20)
