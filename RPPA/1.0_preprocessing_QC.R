#### RPPA analysis: 1.0 pre-processing and QC
#### this script does median centering to normalize the RPPA data plus 
#### doing PCA and correlation analysis to visualize data quality

rm(list = ls())

# Load libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(pheatmap)

# Define input paths
metadata_file <- "samplesheet.csv"
raw_file      <- "RPPA_t.csv" # data is structure as row = raw intensity and column = each sample

# 1. Load and preprocess data

metadata <- read_csv(metadata_file) %>%
  separate(col = sample, sep = "-", into = c("HA", "sample")) %>%
  select(-HA) %>% # clean-up sample names
  mutate(sample = factor(sample, levels = unique(sample)))

raw <- read_csv(raw_file) %>%
  filter(X != "X") %>%
  column_to_rownames("X") %>%
  select(metadata$GMU_ID)

## Median scaling normalization and log2 transform
sample_medians <- apply(raw, 2, median, na.rm = TRUE)
target_median <- median(sample_medians)
scaling_factors <- sample_medians / target_median
raw_scaled <- sweep(raw, 2, scaling_factors, FUN = "/")
raw_scaled_log <- log2(raw_scaled)

## Save pre-processed expression + metadata
saveRDS(list(expr = raw_scaled_log , meta = metadata), file = "results/01_preprocessed_data.rds")

# 2. QC: PCA analysis
pca_data <- t(raw_scaled_log ) %>% as.data.frame()
meta_pca <- metadata %>% column_to_rownames("GMU_ID")

pca <- PCA(pca_data, graph = FALSE)

p_pca <- fviz_pca_ind(pca,
                      geom        = "point",
                      habillage   = meta_pca$sample,
                      addEllipses = FALSE,
                      mean.point  = FALSE,
                      pointshape  = 19,
                      pointsize   = 3,
                      legend.title= "Sample") +
  theme_classic()

print(p_pca)

# 3. QC: Pearson correlation heatmap

corr_mat <- cor(raw_scaled_log, method = "pearson", use = "pairwise.complete.obs")

ann_col <- data.frame(sample = metadata$sample)
rownames(ann_col) <- metadata$GMU_ID        # ensure identical order

cor_plot <-pheatmap(
  corr_mat,
  annotation_col       = ann_col,
  annotation_row       = ann_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method    = "ward.D2",
  display_numbers      = FALSE,
  main                 = "Sample–sample Pearson correlation",
  width                = 6, height = 6
)

print(cor_plot)
