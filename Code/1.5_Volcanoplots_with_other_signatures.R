#### make volcano plot using results from Limma DE analysis and overlap the volcano plot with other known signatures
rm(list = ls())
setwd("directory")

# load necessary packages
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# load ensmbl_gene data base

gene_info_simp <- read.csv("ref/gene_info.csv") %>% mutate(X = NULL) # contain ensembl_gene_id, and gene symbol
gene_info <- read.csv("ref/gene_info_alias.csv") %>% mutate(X = NULL) # ensemble_gene_id, gene symbol, and gene alias

# load RNA-seq data

temp <- read.csv("1_limma/DE_limma_KEAP1_KO_final_241115.csv") %>%
  rename(ensembl_gene_id = X) %>%
  left_join(gene_info_simp, by ="ensembl_gene_id") %>%
  filter(Gene != "" ) 

# filter gene that are significantly regulated (p values < 0.05)
temp.sig <- temp %>%
  filter(adj.P.Val <=0.05)

# load KEAP1-NRF2 gene sets

## PDAC KEAP1-loss signature defined in this study
PKD <- read.csv("3_NRF2up/KNR_signature_240612.csv")  %>%
  mutate(X=NULL) %>%
  left_join(temp, by = "ensembl_gene_id")

## KEAP1 signature published by Singh et al, 2021
SINGH_21 <- read.table("3_NRF2up/signatures/SINGH_2021.txt") %>% 
  left_join(gene_info_simp, by = c("V1" = "Gene")) %>%
  left_join(temp, by = "ensembl_gene_id") %>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))

## KEAP1-NRF2 signature published by Chroley et al, 2012
CHORLEY <- read.csv("3_NRF2up/signatures/Chorley_2012.csv") %>%
  left_join(gene_info_simp, by = "Gene") %>%
  distinct(Gene, .keep_all = T) %>%
  left_join(temp, by = "ensembl_gene_id")  %>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))

## NRF2 signature published by Ibrahim et al, found in Molecular Signature Database
IBRAHIM <- read.csv("3_NRF2up/signatures/IBRAHIM_240617.csv") %>%
  select(-X) %>%
  rename(ensembl_gene_id = ensembl_gene) %>%
  left_join(temp, by = c("ensembl_gene_id"))%>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))

## Curated KEAP1-NRF2 signature published by Morgenstern et al.
Morg_uni <- read.csv("3_NRF2up/signatures/Morgenstern_full.csv") %>%
  left_join(gene_info, by = c("Gene" = "hgnc_symbol")) %>%
  mutate(match_type = ifelse(!is.na(ensembl_gene_id), "hgnc_symbol", NA)) %>%
  {
    unmatched <- filter(., is.na(ensembl_gene_id)) %>%
      select(-any_of(c("ensembl_gene_id", "hgnc_symbol", "external_synonym", "match_type"))) %>%
      left_join(gene_info, by = c("Gene" = "external_synonym")) %>%
      mutate(match_type = ifelse(!is.na(ensembl_gene_id), "external_synonym", NA))
    
    matched <- filter(., !is.na(ensembl_gene_id))
    bind_rows(matched, unmatched)
  } %>%
  select(Gene, ensembl_gene_id) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  left_join(temp, by = "ensembl_gene_id")  %>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))

# KEAP1-NRF2 signature published by Lou et al., 2024
Lou <- read.csv("3_NRF2up/signatures/Lou_2024.csv") %>%
  left_join(gene_info_simp, by = "Gene") %>%
  distinct(Gene, .keep_all = T) %>%
  left_join(temp, by = "ensembl_gene_id")  %>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))


# KEAP1-NRF2 signature published by Liu et al, 2019
## Problem for Liu Gene sets: if genes show up multiple times, an additional "1" will be added to the second one,
## and "2 will be added to the third one

### clean up Liu et al signature
### For rows labeled “nrf2”, remove a trailing “1” *only if* 
### the same gene (without the 1) already appears in the “nrf1” subset.
### Preserve all genuine gene symbols that legitimately end in “1”.

Liu <- read.csv("3_NRF2up/signatures/Liu_2019.csv")  

Liu_fixed <- Liu %>%  #oroginal data
  filter(NRF.member %in% c("nrf1", "nrf2")) %>% 
  { # create a temporary working space                                   
  # Collect the set of bona‑fide gene names found in the nrf1 subset
  nrf1_genes <- unique(.$Gene[.$NRF.member == "nrf1"])
  
  # In the nrf2 subset, drop an added “1” when it is redundant
  mutate(.,
         Gene = if_else(
           NRF.member == "nrf2" &               # operate only on the nrf2 rows
             str_ends(Gene, "1") &                # gene ends with “1”
             str_remove(Gene, "1$") %in% nrf1_genes,  # base name exists in nrf1
           str_remove(Gene, "1$"),              # → strip the trailing “1”
           Gene                                 # otherwise keep the original
         )
  )
}

Liu_UP <- Liu_fixed  %>%
  filter(NRF.member == "nrf2" & Log2.Fold.Change >0 & adjusted.p.value <0.05) %>%
  left_join(gene_info, by = c("Gene" = "hgnc_symbol")) %>%
  mutate(match_type = ifelse(!is.na(ensembl_gene_id), "hgnc_symbol", NA)) %>%
  {
    unmatched <- filter(., is.na(ensembl_gene_id)) %>%
      select(-any_of(c("ensembl_gene_id", "hgnc_symbol", "external_synonym", "match_type"))) %>%
      left_join(gene_info, by = c("Gene" = "external_synonym")) %>%
      mutate(match_type = ifelse(!is.na(ensembl_gene_id), "external_synonym", NA))
    
    matched <- filter(., !is.na(ensembl_gene_id))
    bind_rows(matched, unmatched)
  } %>%
  select(Gene, ensembl_gene_id) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  left_join(temp, by = "ensembl_gene_id")  %>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))

Liu_DN <- Liu_fixed  %>%
  filter(NRF.member == "nrf2" & Log2.Fold.Change <0 & adjusted.p.value <0.05) %>%
  left_join(gene_info, by = c("Gene" = "hgnc_symbol")) %>%
  mutate(match_type = ifelse(!is.na(ensembl_gene_id), "hgnc_symbol", NA)) %>%
  {
    unmatched <- filter(., is.na(ensembl_gene_id)) %>%
      select(-any_of(c("ensembl_gene_id", "hgnc_symbol", "external_synonym", "match_type"))) %>%
      left_join(gene_info, by = c("Gene" = "external_synonym")) %>%
      mutate(match_type = ifelse(!is.na(ensembl_gene_id), "external_synonym", NA))
    
    matched <- filter(., !is.na(ensembl_gene_id))
    bind_rows(matched, unmatched)
  } %>%
  select(Gene, ensembl_gene_id) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  left_join(temp, by = "ensembl_gene_id")  %>%
  filter(!is.na(ensembl_gene_id) & !is.na(logFC))

### consolidate all gene sets

annot_list_2 <-list(
  'SINGH' = SINGH_21 %>% select(ensembl_gene_id, logFC, adj.P.Val),
  'IBRAHIM' = IBRAHIM  %>% select(ensembl_gene_id, logFC, adj.P.Val),
  'CHORLEY' = CHORLEY  %>% select(ensembl_gene_id, logFC, adj.P.Val),
  'MORGENSTERN' = Morg_uni  %>% select(ensembl_gene_id, logFC, adj.P.Val),
  'LIU_UP' = Liu_UP  %>% select(ensembl_gene_id, logFC, adj.P.Val),
  'LIU_DN' = Liu_DN  %>% select(ensembl_gene_id, logFC, adj.P.Val),
  'LOU' = Lou  %>% select(ensembl_gene_id, logFC, adj.P.Val)
)

anno_2 <- bind_rows(annot_list_2, .id = "dataset") %>%
  mutate(direction = ifelse(adj.P.Val > 0.05, 0,
                            ifelse(logFC > 0, 1, -1)))


KEAP1 <-temp %>%
  filter(ensembl_gene_id == "ENSG00000079999")

# make volcano plots
## the PKD genes are plotted in red
## KEAP1 is plotted in purple
p1 <- ggplot(temp) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
  geom_point(data = PKD, aes(x = logFC, y = -log10(adj.P.Val)), fill = "red", color = "black",alpha = 0.9, size = 3, shape = 21) +
  geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
  xlim(-4,4)+
  ylim(-0.5,20)+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0) +
  labs(y="-log(adj. p value)", x="log2FC")+
  theme_classic()+
  theme(axis.text.y   = element_text(size = 15, hjust = 1),
        axis.text.x   = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        axis.title.x  = element_text(size = 15),
        axis.line.y.right = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "none"
  )

print(p1)

# make the tile plots for each signature
## for genes that are not significantly regulated in my RNA-seq --> grey
## for genes that are significantly up-regulated --> red
## for genes taht are significantly down-regulated --> blue

p2 <- ggplot(anno_2, aes(x = logFC, y = dataset, fill = factor(direction))) +
  geom_tile(height = 0.7, width  = 0.01) +
  scale_fill_manual(
    values = c("1"  = "red",    # logFC > 0
               "-1" = "blue",   # logFC < 0
               "3"  = "black"),
    guide = "none"
  ) + 
  xlim(-4,4)+
  scale_y_discrete(limits = rev(names(annot_list_2))) +
  theme_void() +
  theme(
    axis.text.y   = element_text(size = 11, hjust = 1),
    plot.margin    = margin(t = 0, r = 5, b = 5, l = 5) )

print(p3)


# Combine the main volcano plot (p1) with the annotation strip (p2)
combined_plot <- p1 / p2 + 
  plot_layout(heights = c(2.5, 1)) 

# Display the final combined plot
print(combined_plot)

# Save photos
png(paste0("1_limma/volcano_plots_with_other_sigs_5_pvalue_250421.png"), width = 2000, height = 2250, res = 300)

print(combined_plot_2)

dev.off()
