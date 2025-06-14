#### 6.0 DepMap analysis: use DepMap dependency scores database to determine whether topp 200 KEAP1-NRF2 genes
#### are essential genes in PDAC and lung cancer cell lines

rm(list = ls())
# setwd("Directory")

# load Libraries

library(tidyverse)
library(pheatmap)

# Load gene sets & build signature list

## create weighted logFC including the information of p values and logFC
temp <- read.csv("1.1_limma_new/DE_limma_250530_blockingterm.csv") %>%
  select(-X.1, X) %>%
  mutate(
    logFC=-logFC, 
    logp=-log(adj.P.Val),
    scale.logp = (logp - min(logp))/(max(logp)-min(logp)),
    weighted.logFC = logFC * scale.logp,
  ) %>%
  mutate(rank = rank(weighted.logFC))

# filter top 200 up- and donw-regulated genes
PKD_DN <- temp %>%
  filter(rank <=200) 

PKD_UP <- temp %>%
  filter(rank >= (max(temp$rank)-200) )

sig_list <- list(
  KEAP1_KO_DN = PKD_DN %>% select(ensembl_gene_id),
  KEAP1_KO_UP = PKD_UP %>% select(ensembl_gene_id)
)

sig <- bind_rows(sig_list, .id = "dataset") %>%
  left_join(temp %>% select(ensembl_gene_id, Gene), by = "ensembl_gene_id") %>%
  filter(Gene != "") %>%
  mutate(presence = 1)

gene_list <- unique(sig$Gene)

# Load DepMap data & metadata
# median dependency scores for DepMap pancreas and lungs data (pre-downloaded and processed)

depmap <- readRDS("Ref/DepMap_KRAS-MT_PANCREAS_LUNGS.rds") %>%
  select(cell_line, any_of(gene_list)) %>%
  column_to_rownames("cell_line")

depmap_meta <- readRDS("Ref/DepMap_KRAS-MT_PANCREAS_LUNGS_meta.rds")


# Clamp values ±2 & hierarchical clustering (optional)
### NOTE: This may affect hierarchical clustering!!!

mat <- as.matrix(depmap)
# mat[mat >  2.5] <-  2.5
# mat[mat < -2.5] <- -2.5

# row and column clustering
row_hc <- hclust(dist(mat))
col_hc <- hclust(dist(t(mat)))

row_ord <- row_hc$order
col_ord <- col_hc$order

# reorder matrix
mat_ord <- mat[row_ord, col_ord]

# Build annotation for columns (genes)
gene_annotation <- sig %>%
  filter(Gene %in% colnames(mat_ord)) %>%
  distinct(Gene, dataset) %>%
  column_to_rownames("Gene") %>%
  rename(`KEAP1_KO_Group` = dataset)

# Ensure ordering matches heatmap column order
gene_annotation <- gene_annotation[colnames(mat_ord), , drop = FALSE]

#  Define colors for UP and DN gene sets
ann_colors <- list(
  KEAP1_KO_Group = c(
    KEAP1_KO_UP = "red",
    KEAP1_KO_DN = "blue"
  )
)

# Plot heatmap using pheatmap, with annotation of gene expression trends in KEAP1-KO

pheatmap(
  mat_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = gene_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = FALSE,
  fontsize_col = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "DepMap Essentiality: KEAP1_KO Signatures"
)
