<<<<<<< HEAD
#### 1.2 Run differential expression analysis using Limma; use blocking term to block cell lines
=======
#### Run differential expression analysis using Limma
#### use blocking term to block cell lines
>>>>>>> 54959003c7d3a49a37f71f3343907f09b11a1e63

rm(list=ls()) 

# load necessary packages

library(limma)
library(edgeR)
library(tidyverse)

# load data

## raw data with counts
raw <- read.csv("all_counts_noNA_2.csv") %>%
  mutate(X=NULL) %>%
  column_to_rownames("ensembl_gene_id") 

## metadata with experimental conditions
metadata <- read.csv("1_limma/metadata_limma_2.csv")

## gene reference (ensembl_gene_id and gene symbol from biomart)
gene_info <- read.csv("ref/gene_info.csv")

### This checks order
colnames(raw)==metadata$sample

# creat design matrix for limma

group <- factor(paste(metadata$cell_line, metadata$condition, sep="."))
design <- model.matrix(~ 0 + condition, data = metadata)


# Create DGEList object and normalize
dge <- DGEList(counts = raw)
dge <- calcNormFactors(dge)

# Normalization and Voom Transformation (for RNA-seq data):
## Transform counts using voom = make raw counts to log CPM
v <- voom(dge, design, plot = TRUE)

# make cell line blocking term
block <- factor(metadata$cell_line)
corfit <- duplicateCorrelation(v, design, block = block)

# Fit the linear model

fit <- lmFit(v, design, block = block, correlation =  corfit$consensus.correlation)

## KEAP1KO vs. Control

# Define contrast matrix

contrast.matrix_KEAP1KO <- makeContrasts(
  conditionKEAP1KO - conditioncontrol,
  levels=design)

# Apply Contrasts and Compute Statistics:

fit_KEAP1KO <- contrasts.fit(fit, contrast.matrix_KEAP1KO)
fit_KEAP1KO <- eBayes(fit_KEAP1KO)

# Extract result

results_KEAP1KO <- topTable(fit_KEAP1KO, adjust="BH", number=Inf)

results_KEAP1KO <- results_KEAP1KO %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id")

# comment this after saving the results
write.csv(results_KEAP1KO, "1.1_limma_new/DE_limma_250530_blockingterm.csv")

