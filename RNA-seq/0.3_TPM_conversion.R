#### 0.3 covert raw counts to TPM

rm(list = ls())

# load library 
library(tidyverse)

#  Define inputs 
sample_names <- c("SW1990", "PANC1", "Pa16C", "Pa14C_new")
gene_length <- read.csv("ref/gene_lengths_v44_240328.csv")

##########################
# Run TPM conversion 
##########################

## process one sample at a time
## each sample has four conditions, triplicate for each condition = 12 samples in total
for (cell in sample_names) {
  
  # Load raw counts
  raw <- read.csv(paste0("raw_counts/", cell, ".csv")) %>%
    select(-X) %>%
    rename(ensembl_gene_id = Geneid) %>%
    mutate(ensembl_gene_id = gsub("\\..*$", "", ensembl_gene_id)) %>%
    column_to_rownames("ensembl_gene_id")
  
  # Filter low-expressed genes (mean counts ≥ 20 across 12 samples)
  raw_filtered <- raw %>%
    mutate(avg_counts = rowMeans(select(., 1:12), na.rm = TRUE)) %>%
    filter(avg_counts >= 20) %>%
    select(-avg_counts)
  
  # Add prior count to avoid division by zero
  raw_adjusted <- raw_filtered + 0.25
  
  # Merge with gene length and convert to RPK
  temp <- raw_adjusted %>%
    mutate(ensembl_gene_id = rownames(raw_adjusted)) %>%
    left_join(gene_length, by = "ensembl_gene_id") %>%
    mutate(length_kb = length / 1000) %>%
    
  rpk <- sweep(temp[1:12], 1, temp$length, FUN = "/")
  rownames(rpk) <- rownames(raw_adjusted)
  
  # Calculate scaling factor
  scaling_factors <- colSums(rpk) / 1e6
  
  # Final TPM
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")
  
  # Save TPM
  write.csv(tpm, paste0("0_TPM/", cell, "_TPM_with_prior_filter.csv"))
}
#########################################
# Combine and average TPM across samples
#########################################

rm(list = ls())

# load sample names
sample_names <- c("SW1990", "PANC1", "Pa16C", "Pa14C_new")

# prepare an empty dataframe for storing data
all_df <- data.frame(ensembl_gene_id = character())

for (cell in sample_names) {
  
  # Load TPM and extract relevant columns (this selects DMSO and KEAP1-KO data only DMSO (2:4), KEAP1KO (8:10))
  tpm <- read.csv(paste0("0_TPM/", cell, "_TPM_with_prior_filter.csv")) %>%
    rename(ensembl_gene_id = X) %>%
    select(ensembl_gene_id, 2:4, 8:10)  # adjust the selection when needed
  
  # Merge with previous cell line TPMs
  all_df <- full_join(all_df, tpm, by = "ensembl_gene_id")
}

# Save combined TPM data
write.csv(all_df, "0_TPM/all_individual_TPM_KEAP1KOvControl.csv", row.names = FALSE)