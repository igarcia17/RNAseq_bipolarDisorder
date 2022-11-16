setwd("C:/Users/CBM/Desktop/RNAseq_bipolarDisorder")

suppressPackageStartupMessages({
  library(vsn, quietly = TRUE)
  library(gplots, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(cluster, quietly = TRUE)
  library(pheatmap, quietly = TRUE)
  library(grid, quietly = TRUE)
  library(BiocManager, quietly = TRUE)
  library(DESeq2, quietly = TRUE)
  library(AnnotationDbi, quietly = TRUE)
  library('org.Hs.eg.db', quietly = TRUE, character.only = TRUE)
  library(EnhancedVolcano, quietly = TRUE)
  library(tidyr, quietly= TRUE)
})
require("ggrepel", quietly = TRUE)
rm(repos)

#Parameters

cutoff <- 0.05 #significancy p value adjusted
FC_threshold <- 4 #fold change threshold to consider in graphs

##Load data

#Load variables of each sample
config <- './configfile_pedlabels.txt'
sampleTable <- read.table(config, header=TRUE, 
                          colClasses= c('factor','character','factor','factor', 
                                        'factor','factor'))

#Convert the counts into a DeSeq DataSet object
data <- DESeqDataSetFromHTSeqCount(sampleTable, directory=".", 
                                   design = ~ age + gender + PED + condition)

## Analysis
#Pre-filtering: clean some of the noise in the counts
keep <- rowSums(counts(data)) >= 10
data <- data[keep,]
rm(keep)

#With this filter, the object goes from 61806 elements to 28525 elements

# DESeq: not original function to adjust number of iterations
dds <- estimateSizeFactors(data)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit = 10000)
dds_raw <- counts(dds, normalized=FALSE)
dds_normalized <- counts(dds, normalized=TRUE)

#save file with counts and normalized counts
write.table(dds_raw, file="MyResults_DEG/counts_raw.tsv", quote=FALSE, 
             sep = "\t", col.names=NA)
write.table(dds_normalized, file="MyResults_DEG/counts_normalized.tsv", quote=FALSE, 
             sep = "\t", col.names=NA)

#PCA: blind must be FALSE to take into account batch effect
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vst)
mm <- model.matrix(~condition, colData(vst))
mat <- limma::removeBatchEffect(mat,
                                batch=vst$PED, batch2=vst$gender,
                                batch3=vst$age, design=mm)
assay(vst) <- mat
tiff(filename = "MyResults_DEG/PCA.tiff", units="in", width=5, 
     height=5, res=300)
pca <- plotPCA(vst)
pca + ggtitle("Principal Components Plot") + geom_text_repel(aes
                                                             (label=colnames(vst)), 
                                                             size=1)
invisible(dev.off())

#Distances between samples
distRL <- dist(t(mat))
distMat <- as.matrix(distRL)
hc <- hclust(distRL)
hmcol <- colorRampPalette(c("white", "blue"))(299)
tiff(filename = "MyResults_DEG/distances.tiff", units="in", width=5, height=5, res=200)
heatmap.2(distMat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=rev(hmcol),
          margin=c(10, 6), main="Distances matrix", key.title=NA)
invisible(dev.off())

#Dispersion plot
tiff(filename = "MyResults_DEG/dispersion.tiff", units="in", width=5, height=5, res=300)
plotDispEsts(dds, main="Per-gene dispersion estimates")
invisible(dev.off())

##Pair-wise comparisons
#Get factor levels
levels <- unique(sampleTable$condition)
l1 <- toString(levels[2]) #reference level has to be Unaffected
l2 <- toString(levels[1])

suffix <- paste(l1, l2, sep="_vs_")
#Get results
res <- results(dds, contrast=c("condition", l2, l1))
res$FoldChange <- 2^res$log2FoldChange  
res <- res[colnames(res)[c(1,7,2:6)]] # order columns
    
# MAplot
tiff(filename = "MyResults_DEG/maplot.tiff", units="in", width=5, height=5, res=300)
plotMA(res, alpha= 0.05, main=paste("MA-plot", suffix, sep=" "))
invisible(dev.off())

#Add annotation to results
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(res), column="SYMBOL", 
                     keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Hs.eg.db'), keys=row.names(res),
                    column="GENENAME", keytype="ENSEMBL", 
                    multiVals="first") #to obtain description

res <- cbind(symbol, res)
res$description <- description

write.table (res, file="MyResults_DEG/all_genes.tsv", quote=FALSE, sep="\t", col.names=NA)
sig_pval <- subset(res, res$pvalue < cutoff)
write.table (sig_pval, file="MyResults_DEG/sig_pval.tsv", quote=FALSE, sep="\t", col.names=NA)
rm(sig_pval)
#Get most significant genes according to cut off
significant <- subset(res, res$padj < cutoff)
significant <- significant[order(significant$padj),]
#Discard those genes with unbelievable Fold Change (outliers)
significant <- significant[(significant$log2FoldChange >= -FC_threshold) & 
                             (significant$log2FoldChange <= FC_threshold),]
#PCA 2D
interest_genes <- rownames(significant)
dds_sig <- dds[interest_genes,]
vst_sig <- varianceStabilizingTransformation(dds_sig, blind = FALSE)
mat_sig <- assay(vst_sig)
mm_sig <- model.matrix(~condition, colData(vst_sig))
mat_sig <- limma::removeBatchEffect(mat_sig,
                                batch=vst_sig$PED, batch2=vst_sig$gender,
                                batch3=vst_sig$age, design=mm_sig)
assay(vst_sig) <- mat_sig
tiff(filename = "MyResults_DEG/PCA_sig.tiff", units="in", width=5, 
     height=5, res=300)
pca_sig <- plotPCA(vst_sig)
pca_sig + ggtitle("PCA - only significant genes") + geom_text_repel(aes
                                                             (label=colnames(vst_sig)), 
                                                             size=1.5)
invisible(dev.off())

write.table (significant, file="MyResults_DEG/0.05_sig_padj.tsv", quote=FALSE, sep="\t", col.names=NA)

#Volcano plot
discardNA <- !is.na(res$padj)
res2 <- res[discardNA,]
remove_outliers<- (res2$log2FoldChange >= -FC_threshold) & (res2$log2FoldChange <= FC_threshold) 
res2 <- res2[remove_outliers,]
rm(discardNA)
rm(remove_outliers)

jpeg(filename = "MyResults_DEG/VolcanoPLOT.jpeg", units="in", width=8, height=10, res=300)
EnhancedVolcano(res2, lab = res2$symbol, x = 'log2FoldChange', y = 'pvalue',
                pCutoff = 0.0002, FCcutoff= 0.3, #pCutOff is p value for last significant acc to adjp value
                ylim = c(0, 11), xlim = c(-FC_threshold, FC_threshold), labSize = 3,
                legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, title = '', arrowheads = FALSE,
                subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
invisible(dev.off())

#Heatmap
significant001 <- significant[significant$padj < 0.01,]
#Extract from the normalized counts table the data from the significant genes
subcounts <- subset(dds_normalized, rownames(dds_normalized) %in% rownames(significant001))
lsubcounts <- log2(subcounts+1) #added pseudocount 1

#For the plot gene names
sig_symbol <- as.character(significant001$symbol)
make_italics <- function(x){
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}
conditions <- c(l1, l2)
conds <- subset(sampleTable, sampleTable$condition %in% conditions)
samples <- conds$sample

df <- data.frame(condition=conds$condition)
rownames(df) <- samples
my_colour <- list(df=c(l1="skyblue", l2="orange"))

jpeg(filename = "MyResults_DEG/heatmap.jpeg", units="in", width=8, height=5, res=300)
pheatmap(lsubcounts, scale= 'row', cluster_rows = TRUE,
         cluster_cols = TRUE, legend= TRUE, drop_levels = TRUE, 
         labels_row = make_italics(sig_symbol), 
         main = "Heatmap of genes with adjusted p-value < 0.01",
         annotation_col = df, annotation_colors = my_colour,
         treeheight_row = 30, treeheight_col = 20)
invisible(dev.off())

rm(conditions)
rm(conds)
rm(samples)
rm(df)
rm(my_colour)