#### 0. acquire ensembl id and gene symbols from biomart

library(biomaRt)
library(tidyverse)

# select biomart database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene IDs, symbols, and aliases
gene_data <- getBM(attributes = c("ensembl_gene_id", 
                                  "hgnc_symbol", 
                                  "external_synonym",
                                  "description"), 
                   mart = ensembl)

# save data for future use
write.csv(gene_data, "ref/gene_data_full.csv")
