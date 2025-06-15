#### RPPA analysis: 2.0 Differential expression analysis
#### Due to antibody qualities and signal-to-noise raio, heatmap is generally recommend to visualize the result of RPPA
#### Further experimental validation is needed for RPPA analysis

rm(list = ls())

library(tidyverse)
library(pheatmap)

# load differential expression data
DE_all <- read_csv("RPPA_DE_result.csv") %>%
  rename(antibody = X)

condition <- unique(DE_all$contrast)

# filter condition of interest (coi), in this case, "everything vs. vector" is used as examples
coi_vec <- condition[str_detect(condition, "vs_vector")] 

DE_coi <- DE_all %>%
  filter(contrast %in% coi_vec)  # filter the data from coi

## expand differential expression table for heatmap plotting
DE_coi_expand <- DE_coi %>%
  select(antibody, logFC, contrast) %>%
  pivot_wider(values_from = logFC, names_from = contrast) %>%
  column_to_rownames("antibody") # reorganize the table by including logFC and coi

# prepare heatmap plotting

## OPTIONAL: forcing the data to confine withing +1.5 to -1.5 (for visualization purpose)
## CAUTION: This step may affect clustering
DE_capped <- DE_coi_expand
DE_capped[DE_capped >  1.5] <-  1.5
DE_capped[DE_capped < -1.5] <- -1.5

## define the color scale
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

## set color breaks from -2.5 to 2.5 (must match palette length)
my_breaks <- seq(-1.5, 1.5, length.out = 101)  # length = number of colors + 1

# plot heatmap
DE_heatmap <- pheatmap(
  DE_capped,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  color                    = my_palette,
  breaks                   = my_breaks,
  display_numbers          = FALSE,
  main                     = "Heatmap for everything versus Vector",
  fontsize                 = 10,
  cellwidth  = 40, 
  cellheight = 10
)

print(DE_heatmap)

