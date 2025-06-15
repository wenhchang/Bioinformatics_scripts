# 2.0 Run GSEA and Over-representation analysis (ORA)

rm(list = ls())                             

# load necessary packages
library(tidyverse)
library(msigdbr)
library(fgsea)

# load gene info table (ensembl_gene_id, gene_symbol)
gene_info <- read.csv("ref/gene_info.csv")


################
# GSEA
################

# generate empty df for later storage
gsea_H  <- data.frame()
gsea_C2 <- data.frame()

# get gene sets from Molecular Signature Database

## get Hallmark gene sets
H.ensembl.ls <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  select(gs_name, ensembl_gene) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(ensembl_gene), .groups = "drop") %>% 
  deframe()

## get C2 gene sets
C2.ensembl.ls <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  select(gs_name, ensembl_gene) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(ensembl_gene), .groups = "drop") %>% 
  deframe()

# Run GSEA

sample_name <- c("PANC1", "Pa14C_new", "Pa16C", "SW1990")

for (cell in sample_name) {
  
  raw <- read.csv(
    paste0("1_DESEQ2_LFCS_240611/", cell, "_low_count_filtered_lfc_DMSO_mtvwt.csv")
  ) %>% 
    rename(ensembl_gene_id = X)
  
  FC.vec <- raw %>% 
    select(ensembl_gene_id, log2FoldChange) %>% 
    { set_names(.$log2FoldChange, .$ensembl_gene_id) } %>% 
    sort(decreasing = TRUE)
  
  temp_H <- fgseaSimple(
    pathways  = H.ensembl.ls,
    stat      = FC.vec,
    scoreType = "std",
    nperm     = 2000
  ) %>% 
    as.data.frame() %>% 
    mutate(source = cell)
  
  temp_C2 <- fgseaSimple(
    pathways  = C2.ensembl.ls,
    stat      = FC.vec,
    scoreType = "std",
    nperm     = 10000
  ) %>% 
    as.data.frame() %>% 
    mutate(source = cell)
  
  gsea_H  <- rbind(gsea_H,  temp_H)
  gsea_C2 <- rbind(gsea_C2, temp_C2)
}

######################################
# Over-representation analysis (ORA)
######################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGgraph)
library(ReactomePA)

# load gene_info table containing entrez id (not the same gene_info table as above)
gene_info <- read.csv("ref/gene_info_entrez.csv") %>% mutate(X = NULL)

sample_name <- c("PANC1", "Pa14C_new", "Pa16C", "SW1990")

for (cell in sample_name) {
  
  raw <- read.csv(
    paste0("1_DESEQ2_LFCS_240611/", cell, "_low_count_filtered_lfc_DMSO_mtvwt.csv")
  ) %>% 
    rename(ensembl_gene_id = X)
  
  FC <- raw %>% 
    select(ensembl_gene_id, log2FoldChange)
  
  # threshold is set as genes upregulated for more than 50%
  threshold <- log2(1.5)
  
  topTen <- FC %>% 
    filter(log2FoldChange > threshold) %>% 
    left_join(gene_info, by = "ensembl_gene_id")
  
  topTen_list        <- topTen$ensembl_gene_id
  topTen_list.entrez <- topTen$entrezgene_id
  
  # ORA for GO gene sets
  enrich_go <- enrichGO(
    gene          = topTen_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  ) %>% 
    as.data.frame() %>% 
    rename(category = ONTOLOGY) %>% 
    mutate(db = "GO") %>% 
    select(db, category, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)
  
  # ORA for KEGG gene sets
  enrich_KEGG <- enrichKEGG(
    gene          = topTen_list.entrez,
    organism      = "hsa",
    keyType       = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  ) %>% 
    as.data.frame() %>% 
    mutate(db = "KEGG") %>% 
    select(db, category, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)
  
  # ORA for Reactome gene sets
  enrich_Reactome <- enrichPathway(
    gene          = topTen_list.entrez,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  ) %>% 
    as.data.frame() %>% 
    mutate(db = "REACTOME", category = "") %>% 
    select(db, category, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)
  
  # combine all results
  Total_ORA <- rbind(enrich_go, enrich_KEGG, enrich_Reactome) %>%
    mutate(sample = cell)
  
  write.csv(
    Total_ORA,
    paste0("2_GSEA_240611/ORA_", cell, "_FC1.5_240621.csv")
  )
}