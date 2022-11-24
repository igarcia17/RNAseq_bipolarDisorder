#Weighted Gene Coexpression Network Analysis part 4, for identification
#of driver genes and comparison with previous DGE results.

suppressPackageStartupMessages({
  library(rstudioapi, quietly=T)
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
  library(VennDiagram, quietly = T)
})
options(stringsAsFactors = F)
enableWGCNAThreads()

#Paths to files
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

filesD <- 'results_WGCNA/'
inputData <- paste0(filesD,'initialData_datExpr_sampleTable.RData')
inputNet <- paste0(filesD, 'network_manual_construction.RData')
inputSigMods <- paste0(filesD, 'summary_sigMods.RData')
inputSigGenes <- 'results_DGE/0.05_sig_padj.tsv'
load(inputData)
load(inputNet)
load(inputSigMods)

vennF <- paste0(filesD, 'common_modules_fromWGCNA_DEGs.tiff')
pieF <- paste0(filesD, 'pie_chart_common_modules.tiff')

#Relate to previous DEG results
resDGE <- read.table(file = inputSigGenes, sep = '\t', header = TRUE, row.names=1)
resDGE <- resDGE[order(row.names(resDGE)),]
DEGs <- rownames(resDGE)

DEGinfo <- filter(geneInfo, rownames(geneInfo) %in% DEGs)
DEGinfo <- DEGinfo[order(row.names(DEGinfo)),]

#Why are there 4 genes missing?
missingGenes <- setdiff(DEGs, rownames(DEGinfo))
missGenesDGEs <- filter(resDGE, rownames(resDGE) %in% missingGenes)
#Looking at baseMean, we can see that they have probably been filtered out for the 
#WGCNA because of low counts, or because they are not significant enough

#What significant modules appear in the DEGs? Venn diagram
DEGmods_list <- unique(DEGinfo$moduleColor)
x <- list('Modules of significant genes' = DEGmods_list, 
          'Significant modules' = sigMods)
vennPlot <- venn.diagram(x, disable.logging = T, print.mode = 'percent',
                  fill = c("orange", "blue"), filename= NULL,
                  alpha = c(0.5, 0.5), cat.cex = 0.9, cex=1.3,cat.pos = 3)
vennPlot[[5]]$label  <- paste(setdiff(DEGmods_list, sigMods), collapse="\n")
vennPlot[[6]]$label <- paste(setdiff(sigMods, DEGmods_list)  , collapse="\n")
vennPlot[[7]]$label <- paste(intersect(DEGmods_list, sigMods), collapse="\n")  

tiff(file = vennF)
grid.newpage()
grid.draw(vennPlot)
invisible(dev.off())

#What percentage of DEGS are in a significant module?
DEGmods_df <- as.data.frame(DEGinfo$moduleColor)
colnames(DEGmods_df) <- 'WGCNA_module'
noWGCNAsig <- !(DEGmods_df$WGCNA_module %in% sigMods)
DEGmods_df$WGCNA_module[noWGCNAsig] <- 'Other'

DEGmods_summary <- DEGmods_df %>% 
  group_by(WGCNA_module) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

tiff(filename = pieF)
ggplot(DEGmods_summary, aes(x="", y=perc, fill=WGCNA_module)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=180) +
  theme_void() +
  geom_text(aes(x = 1.35, label = labels),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
invisible(dev.off())
