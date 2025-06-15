#### 4.0 Use Euler Diagram to compare the overlap of KEAP1-NRF2 gene signature with 
#### KRAS-, MYC-, or YAP/TEAD-dependent signature

rm(list=ls()) # remove all unnecessary intermediate data
setwd("D:/Dropbox/Der lab/RNA seq/All_data_recalculation")

# load necessary package
library(tidyverse)
library(msigdbr)

# load reference table (ensembl_gene_id, alias, and gene symbol)

gene_info <- read.csv("ref/gene_info.csv") %>%
  select(-X)

gene_data <- read.csv("ref/gene_data.csv")

# load signature tables
## KRAS-dependent: Klomp et al, Science 2024
## MYC-dependent: Hallmark MYC V1 and V2 combined
## YAP/TEAD-dependent: Pham et al, Cancer Discovery 2021 (Cluster 2)

## load PKD signature (PDAC KEAP1-loss signautre)
PKD <- read.csv("3_NRF2up/KNR_signature_240612.csv") %>% 
  select (-X) %>%
  left_join(gene_info, by = "ensembl_gene_id")
PKD.vec <- PKD$ensembl_gene_id

# load Hallmark MYC signature
H <-msigdbr(species = "Homo sapiens", category = "H") #load gene sets
H.ensembl <- H %>%
  select(gs_name, ensembl_gene) 
H_MYC <- H.ensembl %>%
  filter(grepl("MYC", gs_name)) %>%
  distinct(ensembl_gene, .keep_all = T)

## pull all Hallmark MYC genes together, remove overlap
## check overlap make sure all genes in MYC V1 and V2 are unique
# H_MYC_V1 <- H.ensembl %>%
#   filter(gs_name == "HALLMARK_MYC_TARGETS_V1")
# 
# H_MYC_V2 <- H.ensembl %>%
#   filter(gs_name == "HALLMARK_MYC_TARGETS_V2")
# 
# overlap <- intersect(H_MYC_V1$ensembl_gene, H_MYC_V2$ensembl_gene)

H_MYC.vec <- H_MYC$ensembl_gene

# KRAS signature (Klomp et.al, 2024)

KRAS <- read.csv("other_data/KRAS_up.csv") %>%
  mutate(ensembl_gene_id = gsub("\\..*$", "", ensembl_gene_id))
KRAS.vec <- KRAS$ensembl_gene_id

# Pham's YAP/TEAD signature (Pham et al, 2021)

YAP <- read.csv("other_data/Pham_cluster2.csv") %>%
  select(-term) %>%
  left_join(gene_info, by= c("gene" = "Gene")) %>%
  distinct(gene, .keep_all = T) # Pham signature cluster 2

YAP_unmatch <- YAP %>%
  filter(is.na(ensembl_gene_id)) %>%
  select(gene) %>%
  mutate(gene = toupper(gene)) %>%
  left_join(gene_data, by = c("gene"="external_synonym")) %>%
  select(gene, ensembl_gene_id) %>%
  filter(!is.na(ensembl_gene_id))
# make sure all gene symbol are corresponded to their ensembl gene id

YAP <- YAP %>%
  left_join(YAP_unmatch, by = "gene") %>%
  mutate(ensembl_gene_id = coalesce(ensembl_gene_id.x, ensembl_gene_id.y)) %>%
  filter(!is.na(ensembl_gene_id))

YAP.vec <- YAP$ensembl_gene_id

# Prepare Euler diagram

library(eulerr)
# Create a named list of your sets
Euler_data <- list(
  PKD = PKD.vec,
  MYC = H_MYC.vec,
  KRAS = KRAS.vec,
  YAP = YAP.vec
)
# Professional color palette (colorblind-friendly)
colors <- c("purple3", "steelblue","cadetblue3","violetred")

# Plot the Euler diagram with professional colors
fit <- euler(Euler_data)
p1 <- plot(fit, 
           fills = list(fill = colors, alpha = 0.7),  # Apply colors and set transparency
           edges = list(col = "black", lty = 1.5, lwd = 1.5),  # Add black edges for clarity
           labels = list(col = "black", fontsize = 18))  # Customize label appearance

print(p1)

png("10.5_Euler/HMYC_KRAS_KEAP1_YAP_signature_250327.png", width = 2000, height = 2000, res = 600)

print(p1)

dev.off()


# plot Euler diagram (without lable)

p2 <- plot(fit,
           fills = list(fill = colors, alpha = 0.4),  # Apply colors and set transparency
           edges = list(col = "black", lty = 1, lwd = 1),  # Add black edges for clarity
           labels = FALSE) #

png("10.5_Euler/HMYC_KRAS_KEAP1_YAP_signature_NOLABEL_250327.png", width = 2000, height = 2000, res = 600)

print(p2)

dev.off()

overlap <- fit$original.values %>% as.data.frame()

# save Euler analysis for future use
write.csv(overlap, "10.5_Euler/Euler_HMYC_KRAS_KEAP1_YAP.csv")


