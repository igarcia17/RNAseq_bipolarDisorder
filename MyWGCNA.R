setwd("C:/Users/CBM/Desktop/WGNCA/counts")

library(WGCNA)
options(stringsAsFactors = FALSE)

#1. Load data
df <- read.table(file = 'MyResults_DEG/normalized_counts.tsv', sep = '\t', header = FALSE)

config <- './configfile.txt'
sampleTable <- read.table(config, header=TRUE, colClasses= c('factor','character',
                                                             'factor',
                                                             'factor', 'factor',
                                                             'factor'))
#2. Network construction and module detection  

#3. Relating modules to external clinical traits and identifying important genes

#4. Interfacing network analysis with other data such as functional annotation and gene ontology

#5. Network visualization using WGCNA functions
