#Run a Differentially Expressed Genes Analysis usng as input the count files 
#produced in previous step

suppressPackageStartupMessages({
  library(rstudioapi)
  library(gplots, quietly = T)
  library(ggplot2, quietly = TRUE)
  library(pheatmap, quietly = TRUE)
  library(DESeq2, quietly = TRUE)
  library('org.Hs.eg.db', quietly = TRUE, character.only = TRUE)
  library(EnhancedVolcano, quietly = TRUE)
})
require("ggrepel", quietly = TRUE)

#Paths
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
#Input
configFile <- 'configfile.txt'
#Outputs
resD <- 'results_DGE/'
rawCountsF <- paste0(resD,"counts_raw.tsv")
normCountsF <- paste0(resD,"counts_normalized.tsv")
PCAF <- paste0(resD,"PCA.jpeg")
distancesF <- paste0(resD,"distances.jpeg")
dispersionF <- paste0(resD,"dispersion.tiff")
MAplotF <- paste0(resD,"maplot.jpeg")
genesTSV <- paste0(resD,"all_genes.tsv")
sigTSV <- paste0(resD,"sig_pval.tsv")
sigPCAF <- paste0(resD,"PCA_sig.jpeg")
alphasigTSV <- paste0(resD,"0.05_sig_padj.tsv")
volcanoF <- paste0(resD,"volcanoPlot.jpeg")
heatmapF <- paste0(resD,"heatmap.jpeg")
DESEqResultsF <- paste0(resD, 'deseq_objects.RData')

#Parameters
cutoff <- 0.05 #significancy p value adjusted
FCthres <- 4 #fold change threshold to consider in graphs
covs <- T

#Functions
make_italics <- function(x){
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

#Load variables of each sample
sampleTable <- read.table(configFile, header=TRUE, 
                          colClasses= c('factor','character','factor','factor', 
                                        'factor','factor'))

#Convert the counts into a DeSeq DataSet object
data <- DESeqDataSetFromHTSeqCount(sampleTable, directory=".", 
                                   design = ~ age + gender + PED + condition)

## Analysis
#Pre-filtering: clean some of the noise in the counts
keep <- rowSums(counts(data)) >= 10
data <- data[keep,]
#With this filter, the object goes from 61806 elements to 28525 elements

# DESeq: original DESEQ() function doesnt allow to adjust number of iterations
dds <- estimateSizeFactors(data)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit = 10000)
dds_raw <- counts(dds, normalized=FALSE)
dds_normalized <- counts(dds, normalized=TRUE)

#save file with counts and normalized counts
write.table(dds_raw, file=rawCountsF, quote=FALSE, 
             sep = "\t", col.names=NA)
write.table(dds_normalized, file=normCountsF, quote=FALSE, 
             sep = "\t", col.names=NA)
#Save DESeq results for GSEA
save(dds, file = DESEqResultsF)

#PCA: blind must be FALSE to take into account batch effect
if (covs) {
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)
  mat <- assay(vst)
  mm <- model.matrix(~condition, colData(vst))
  mat <- limma::removeBatchEffect(mat,
                                batch=vst$PED, batch2=vst$gender,
                                batch3=vst$age, design=mm)
  assay(vst) <- mat
} else {
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)
  mat <- assay(vst)
}

jpeg(filename = PCAF, width=900, height=900, quality=300)
pca <- plotPCA(vst)
title <- "Principal Components Plot"
pca + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst)), size=5, point.padding = 0.6)
invisible(dev.off())

#Distances between samples
distRL <- dist(t(mat))
distMat <- as.matrix(distRL)
hc <- hclust(distRL)
hmcol <- colorRampPalette(c("white", "blue"))(299)
title <- "Distances matrix"
jpeg(filename = distancesF, width=900, height=900, quality=300)
heatmap.2(distMat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=rev(hmcol),
          margin=c(10, 6), main=title, key.title=NA)
invisible(dev.off())

#Dispersion plot
tiff(filename = dispersionF, units="in", width=5, height=5, res=300)
title <- "Per-gene dispersion estimates"
plotDispEsts(dds, main=title)
invisible(dev.off())

##Pair-wise comparisons
#Get factor levels
levels <- unique(sampleTable$condition)
l1 <- toString(levels[2]) #reference level has to be Unaffected
l2 <- toString(levels[1])
suffix <- paste(l1, l2, sep="_vs_")

#Get results
res <- results(dds, contrast=c("condition", l2, l1))
res$FoldChange <- 2^res$log2FoldChange  #have actual fold change
res <- res[colnames(res)[c(1,7,2:6)]] # order columns
    
# MAplot
tiff(filename = MAplotF, units="in", width=5, height=5, res=300)
title <- paste("MA-plot", suffix, sep=" ")
plotMA(res, alpha= 0.05, main=title, colSig = 'red')
invisible(dev.off())


jpeg(filename = MAplotF, width=900, height=900, quality=300)
title <- paste("MA-plot", suffix, sep=" ")
plotMA(res, alpha= cutoff, main=title, colSig = 'red', cex = 1.2)
invisible(dev.off())

#Add annotation to results
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(res), column="SYMBOL", 
                     keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Hs.eg.db'), keys=row.names(res),
                    column="GENENAME", keytype="ENSEMBL", 
                    multiVals="first") #to obtain description

res <- cbind(symbol, res)
res$description <- description

write.table (res, file=genesTSV, quote=FALSE, sep="\t", col.names=NA)
sig_pval <- subset(res, res$pvalue < cutoff)
write.table (sig_pval, file=sigTSV, quote=FALSE, sep="\t", col.names=NA)

#Get most significant genes according to cut off
significant <- subset(res, res$padj < cutoff)
significant <- significant[order(significant$padj),]
#Discard those genes with unbelievable Fold Change (outliers)
significant <- significant[(significant$log2FoldChange >= -FCthres) & 
                             (significant$log2FoldChange <= FCthres),]

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

jpeg(filename = sigPCAF, width=900, height=900, quality=300)
pca_sig <- plotPCA(vst_sig)
title <- "PCA - only significant genes"
pca_sig + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst_sig)), size=5, point.padding = 0.6)
invisible(dev.off())

write.table (significant, file=alphasigTSV, quote=FALSE, sep="\t", col.names=NA)

#Volcano plot
discardNA <- !is.na(res$padj)
res2 <- res[discardNA,]
remove_outliers<- (res2$log2FoldChange >= -FCthres) & (res2$log2FoldChange <= FCthres) 
res2 <- res2[remove_outliers,]

jpeg(filename = volcanoF, units="in", width=8, height=10, res=300)
EnhancedVolcano(res2, lab = res2$symbol, x = 'log2FoldChange', y = 'pvalue',
                pCutoff = 0.0002, FCcutoff= 0.3, 
                #pCutOff is p value for last significant acc to adjp value
                ylim = c(0, 11), xlim = c(-FCthres, FCthres), labSize = 3,
                legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, title = '', arrowheads = FALSE,
                subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
invisible(dev.off())

#Heatmap
significant001 <- significant[significant$padj < 0.01,]
#Extract from the normalized counts table the data from the significant genes
subcounts <- subset(dds_normalized, rownames(dds_normalized) %in% 
                      rownames(significant001))
lsubcounts <- log2(subcounts+1) #added pseudocount 1

#For the plot gene names
sig_symbol <- as.character(significant001$symbol)
conditions <- c(l1, l2)
conds <- subset(sampleTable, sampleTable$condition %in% conditions)
samples <- conds$sample

df <- data.frame(condition=conds$condition)
rownames(df) <- samples
my_colour <- list(df=c(l1="skyblue", l2="orange"))
title <- "Heatmap of genes with adjusted p-value < 0.01"

jpeg(filename = heatmapF, units="in", width=8, height=5, res=300)
pheatmap(lsubcounts, scale= 'row', cluster_rows = TRUE,
         cluster_cols = TRUE, legend= TRUE, drop_levels = TRUE, 
         labels_row = make_italics(sig_symbol), 
         main = title,
         annotation_col = df, annotation_colors = my_colour,
         treeheight_row = 30, treeheight_col = 20)
invisible(dev.off())
