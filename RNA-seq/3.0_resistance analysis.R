#### 3.0 Resistance analysis

##############################################################################
#  Determine whether activation of KEAP1/NRF2 correlates with 
#  resistance to KRAS inhibitors
#  Methods: 
#  1. determine the single cell gene expression scores using "Singscore"
#  2. plot the expression scores vs. drug response to see if the resistant
#     models tend to have higher expression scores
#
#  Sections                                                            
#    1. Adagrasib-treated patients (clinical samples)                  
#    2. Sotorasib-treated PDX                                          
#    3. CCLE G12D mutant cell lines                                    
#    4. Julien et al. resistant vs parental cell lines                 
#############################################################################

# Load libraries
library(tidyverse)   # ggplot2 / dplyr / tidyr / tibble / readr
library(singscore)   # rankGenes(), simpleScore()
library(fgsea)       # fgseaSimple() – used in section 4

# load reference data
gene_info <- read.csv("../ref/gene_info.csv") %>%
  mutate(X=NULL)

gene_length <- read.csv("../Ref/gene_lengths_v44_240328.csv") 


# Load other KEAP1-NRF2 signatures

Morg_uni <- read.csv("../3_NRF2up/signatures/Morg_uni.csv") %>%
  select(-X) %>%
  mutate(gs_name = 'MORGENSTERN') 
Morg_uni.list <- Morg_uni$ensembl_gene_id

Liu_UP <- read.csv("../3_NRF2up/signatures/Liu_UP.csv")%>%
  select(-X) %>%
  mutate(gs_name = 'LIU_UP') 
Liu_UP.list <- Liu_UP$ensembl_gene_id 

Liu_DN <- read.csv("../3_NRF2up/signatures/Liu_DN.csv")%>%
  select(-X) %>%
  mutate(gs_name = 'LIU_DN') 
Liu_DN.list <- Liu_DN$ensembl_gene_id

Lou <- read.csv("../3_NRF2up/signatures/Lou_2024.csv") %>%
  left_join(gene_info, by = "Gene") %>%
  distinct(Gene, .keep_all = T) 
Lou.list <- Lou$ensembl_gene_id

IBRAHIM <- read.csv("../3_NRF2up/signatures/IBRAHIM_240617.csv") %>% 
  select(ensembl_gene) %>%
  mutate(X=NULL, gs_name = 'IBRAHIM') %>%
  rename(ensembl_gene_id = ensembl_gene) %>%
  select(gs_name, everything()) 

IBRAHIM_list <- IBRAHIM$ensembl_gene_id

Chorley <- read.csv("../3_NRF2up/signatures/Chorley_2012.csv") %>% 
  left_join(gene_info, by ="Gene")%>%
  select(ensembl_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = T)%>%
  mutate(X=NULL, gs_name = "Chorley")%>%
  select(gs_name, everything())
Chorley_list <- Chorley$ensembl_gene_id

AS_2021 <- read.csv("../3_NRF2up/signatures/AS_2021.csv") %>% 
  left_join(gene_info, by ="Gene")%>%
  select(ensembl_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = T)%>%
  mutate(X=NULL, gs_name = "SINGH")%>%
  select(gs_name, everything())

AS_2021_list <- AS_2021$ensembl_gene_id

# old version of PKD signature
PKD.old <- read.csv("../3_NRF2up/KNR_signature_240612.csv") %>% 
  mutate(X=NULL, gs_name = "PKD.old") 
PKD.old.ls <- PKD.old$ensembl_gene_id

# new version of PKD signature

PKD.new <- read.csv("NEW_SIGNATURE") %>% 
  mutate(X=NULL, gs_name = "PKD.new") 
PKD.new.ls <- PKD.new$ensembl_gene_id

# direct NRF2 targets signature

direct.NRF2 <- read.csv("dir_NRF2_SIGNATURE") %>% 
  mutate(X=NULL, gs_name = "direct_NRF2") 
direct.NRF2.ls <- direct.NRF2$ensembl_gene_id

#######################################################################
# 1. Adagrasib-treated patients
### Compute singscore for each published/new signature
### Plot screening-timepoint samples coloured by response
#######################################################################

# Import expression matrix (rows = genes, cols = samples) 
G12Ci_patient <- read.csv("../7_Singscore/G12Ci_patient_data_240329.csv") %>%
  column_to_rownames("X")

# Prepare gene-rank matrix 
ranked_data <- rankGenes(G12Ci_patient)

# Run singscore for each signature 
list_of_vectors <- list(
  Chorley      = Chorley_list,
  Singh        = AS_2021_list,
  Ibrahim      = IBRAHIM_list,
  PKD.new      = PKD.new.ls,
  direct_NRF2  = direct.NRF2.ls,
  PKD.old      = PKD.old.ls
)

exp_score_result <- data.frame()

for (sig_name in names(list_of_vectors)) {
  sig_genes <- list_of_vectors[[sig_name]]
  
  tmp <- simpleScore(ranked_data, upSet = sig_genes) %>%
    mutate(gs_name   = sig_name,
           sample_ID = rownames(.))
  
  exp_score_result <- rbind(exp_score_result, tmp)
}

# Reshape & annotate response status 
exp_score_wide <- exp_score_result %>%
  select(sample_ID, TotalScore, gs_name) %>%
  pivot_wider(names_from = gs_name, values_from = TotalScore) %>%
  mutate(suffix    = sub(".*\\.", "", sample_ID),
         sample_ID = sub("\\.[^.]*$", "", sample_ID))

screening_only <- exp_score_wide %>%
  filter(suffix == "Screening")

plot_order <- c("PKD.new", "direct_NRF2", "PKD.old",
                "Singh", "Chorley", "Ibrahim")

data_plot <- screening_only %>%
  pivot_longer(cols = Chorley:PKD.old, values_to = "exp_score") %>%
  mutate(name     = factor(name, levels = plot_order),
         response = ifelse(sample_ID %in% c("X001.806.018", "X001.814.001"),
                           "Sensitive", "Resistant"))

# Dot-plot 
p_adagrasib <- ggplot(data_plot) +
  geom_point(aes(x = exp_score, y = name, fill = response),
             shape = 21, stroke = 1, size = 5) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c(Resistant = "blue", Sensitive = "red")) +
  labs(x = "Singscore", y = NULL,
       title = "Adagrasib – screening samples") +
  theme(axis.text.y  = element_text(size = 24),
        axis.text.x  = element_text(size = 17, angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10),
        plot.background   = element_rect(fill = NA, colour = NA),
        panel.background  = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA, colour = NA))

print(p_adagrasib)



#######################################################################
# 2. Sotorasib-treated PDX
#######################################################################

# Import counts & convert to TPM 
G12Ci_PDX <- read.table("../7_Singscore/GSE204753_counts.txt", header = TRUE) %>%
  left_join(gene_info,  by = c("Name" = "Gene"))               %>%
  distinct(Name, .keep_all = TRUE)                              %>%
  left_join(gene_length, by = "ensembl_gene_id")                %>%
  filter(!is.na(length))                                        %>%
  mutate(length = length / 1000)                                %>%
  column_to_rownames("ensembl_gene_id")                          %>%
  select(-Name)

## RPK → TPM
rpk      <- sweep(G12Ci_PDX[1:9], 1, G12Ci_PDX$length, FUN = "/")
sum_RPK  <- colSums(rpk) / 1e6
TPM      <- sweep(rpk,  2, sum_RPK, FUN = "/")

# Run Singscore analysis
ranked_data <- rankGenes(TPM)

exp_score_result <- data.frame()

for (sig_name in names(list_of_vectors)) {
  sig_genes <- list_of_vectors[[sig_name]]
  
  tmp <- simpleScore(ranked_data, upSet = sig_genes) %>%
    mutate(gs_name   = sig_name,
           sample_ID = rownames(.))
  
  exp_score_result <- rbind(exp_score_result, tmp)
}

exp_score_wide <- exp_score_result %>%
  select(sample_ID, TotalScore, gs_name) %>%
  pivot_wider(names_from = gs_name, values_from = TotalScore) %>%
  filter(sample_ID != "length")

plot_order <- c("PKD.new", "direct_NRF2", "PKD.old",
                "Singh", "Chorley", "Ibrahim")

data_plot <- exp_score_wide %>%
  pivot_longer(cols = Chorley:PKD.old, values_to = "exp_score") %>%
  mutate(name     = factor(name, levels = plot_order),
         response = ifelse(sample_ID %in% c("CTL1", "CTL2", "CTL3"),
                           "Sensitive", "Resistant"))

# Dot-plot 
p_PDX <- ggplot(data_plot) +
  geom_point(aes(x = exp_score, y = name, fill = response),
             shape = 21, stroke = 1, size = 5) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c(Resistant = "blue", Sensitive = "red")) +
  labs(x = "Singscore", y = NULL,
       title = "Sotorasib PDX") +
  theme(axis.text.y  = element_text(size = 24),
        axis.text.x  = element_text(size = 17, angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10),
        plot.background   = element_rect(fill = NA, colour = NA),
        panel.background  = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA, colour = NA))

print(p_PDX)



#######################################################################
# 3. CCLE G12D mutant cell lines
#    – Export singscore + MRTX1133 response for plotting in Prism
#######################################################################

# Import TPM & rank
G12Di_TPM <- read.csv("../Julien_RNAseq/Cell_line/CCLE_G12Di_TPM.csv") %>%
  select(-X) %>%
  column_to_rownames("ensembl_gene_id") 
## data was pre-downloaded from CCLE database, filter only KRAS(G12D)-mutant cancer

ranked_data <- rankGenes(G12Di_TPM)

# Run Singscore 
## use direct NRF2 targets (change to other signature if needed)
result <- simpleScore(ranked_data, upSet = direct_NRF2.ls) %>%
  rownames_to_column("CCLE_name")

## 3C  Merge with drug-response metadata --------------------------------
MRTX1133_res <- read.csv("../Julien_RNAseq/Cell_line/MRTX1133_response.csv") %>%
  select(CCLE_name, G12D_response = G12D_resoponse) %>%
  filter(G12D_response != "") %>%
  left_join(result, by = "CCLE_name")

# write.csv(MRTX1133_res,"MRTX1133_res_direct_NRF2.csv", row.names = FALSE)
## this results are plotted separately using Prism GraphPad


#######################################################################
# 4. Julien’s resistant vs parental cell lines
#######################################################################

# Import DE data of PANC1 and Panc0203, MRTX1133-resistant cells vs. parental counterpart
# data was published by Dilly et al, 2024

## bind the data together for future analysis
all.data <- bind_rows(
  read.csv("../Julien_RNAseq/GSE269985_RAW/Human_cell_pathway_analysis/PANC1_resistance.csv") %>%
    rename(ensembl_gene_id = X) %>%
    mutate(cell = "PANC1"),
  read.csv("../Julien_RNAseq/GSE269985_RAW/Human_cell_pathway_analysis/Panc0203_resistance.csv") %>%
    rename(ensembl_gene_id = X) %>%
    mutate(cell = "Panc0203")
)

# Assemble pathway list 
all.gs <- bind_rows(
  # Hallmark,   # uncomment to include Hallmark
  IBRAHIM     %>% select(gs_name, ensembl_gene_id),
  AS_2021     %>% select(gs_name, ensembl_gene_id),
  Chorley     %>% select(gs_name, ensembl_gene_id),
  PKD.new  %>% select(gs_name = Signature, ensembl_gene_id),
  direct.NRF2 %>% select(gs_name, ensembl_gene_id),
  PKD.old     %>% select(gs_name, ensembl_gene_id)
)

all.gs_ls <- all.gs %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(ensembl_gene_id), .groups = "drop") %>%
  deframe()

# Run GSEA for each cell line 
gsea_result <- data.frame()

for (cell_id in c("PANC1", "Panc0203")) {
  
  FC <- all.data %>%
    filter(cell == cell_id) %>%
    select(ensembl_gene_id, log2FoldChange) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  FC_vec <- set_names(FC$log2FoldChange, FC$ensembl_gene_id)
  
  tmp <- fgseaSimple(
    pathways  = all.gs_ls,
    stat      = FC_vec,
    scoreType = "std",
    nperm     = 2000
  ) %>%
    mutate(source = cell_id)
  
  gsea_result <- rbind(gsea_result, tmp)
}

# Format & plot 
gsea_result <- gsea_result %>%
  select(pathway, padj, NES, source) %>%
  filter(!str_detect(pathway, "HALLMARK_")) %>%
  mutate(pathway = recode(pathway,
                          "PKD_UP.new" = "PKD.new",
                          "SINGH"      = "Singh",
                          "IBRAHIM"    = "Ibrahim"))

plot_order <- c("Ibrahim", "Chorley", "Singh",
                "PKD.old", "PKD.new", "direct_NRF2")

p_res_vs_par <- ggplot(gsea_result,
                       aes(y = pathway, x = source,
                           fill = NES, size = -log10(padj),
                           color = ifelse(padj > 0.05, "gray75", "black"))) +
  geom_point(shape = 21, stroke = 1) +
  scale_fill_gradient(low = "white", high = "red4",
                      limits = c(0, 2), name = "NES") +
  scale_y_discrete(limits = plot_order) +
  scale_size(range = c(0, 10)) +
  scale_color_identity() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title  = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    plot.background   = element_rect(fill = NA, colour = NA),
    panel.background  = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA)
  )

print(p_res_vs_par)