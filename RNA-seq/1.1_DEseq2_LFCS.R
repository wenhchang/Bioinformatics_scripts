#### 1.1 Run differential expression analysis using DESeq2 with log2 fold change shrinkage

rm(list = ls())

# Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)

# Parameters 
cell <- "Pa14C_new"
count_path <- paste0("raw_counts/", cell, "_filter_low_counts_240611.csv")
meta_path  <- paste0("metadata/", cell, "_meta_simple.csv")
gene_info_path <- "gene_info.csv"
output_dir <- "results/01_DESEQ2_LFCS_240611/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data 
gene_info <- read.csv(gene_info_path)

rawcounts <- read.csv(count_path) %>%
  column_to_rownames("X")

metadata <- read.csv(meta_path) %>%
  select(-X)
rownames(metadata) <- as.character(unlist(metadata[,1]))

# Check sample order consistency 
all(rownames(metadata) == colnames(rawcounts))

# Set reference level for condition 
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "wt_DMSO")

# Construct DESeq2 dataset 
dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = metadata,
                              design = ~condition) %>%
  estimateSizeFactors()

# QC: Sample correlation heatmap 
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)

cor_plot <- pheatmap(vsd_cor, annotation_col = metadata["condition"],
         main = paste(cell, "- Sample Correlation"))

print(cor_plot)

# QC: PCA analysis
pca_plot <- plotPCA(vsd, intgroup = "condition")
print(pca_plot)

# Run DE analysis 
dds_DE <- DESeq(dds)

# Extract results: mt_DMSO (KEAP1-KO + DMSO) vs wt_DMSO (WT + DMSO)
res <- results(dds_DE, 
               contrast = c("condition", "mt_DMSO", "wt_DMSO"),
               alpha = 0.05)

# Perform log2 fold change shrinkage using 'ashr' method
## This step improves the estimation of log2FC for low-count genes,
## reducing variability and making downstream visualization (e.g. volcano plots) more stable.

res_lfc <- lfcShrink(dds_DE,
                     contrast = c("condition", "mt_DMSO", "wt_DMSO"),
                     type = "ashr",
                     res = res) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id")

# Save result 
output_file <- paste0(output_dir, cell, "_lfc_DMSO_mtvwt.csv") # change the name based on the condition
write.csv(res_lfc, output_file, row.names = FALSE)