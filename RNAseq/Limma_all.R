#### 1.1 limma
#### goal: to calculate the most accurate KEAP1-NRF2 signature
#### use blocking term to block cell lines

rm(list=ls()) # remove all unnecessary intermediate data
# setwd("D:/Dropbox/Der lab/RNA seq/All_data_recalculation")
setwd("~/Library/CloudStorage/Dropbox/Der lab/RNA seq/All_data_recalculation")

#load library
library(limma)
library(edgeR)
library(tidyverse)

# load data

raw <- read.csv("all_counts_noNA_2.csv") %>%
  mutate(X=NULL) %>%
  column_to_rownames("ensembl_gene_id") 

metadata <- read.csv("1_limma/metadata_limma_2.csv")

gene_info <- read.csv("ref/gene_info.csv")

# This checks order
colnames(raw)==metadata$sample

# creat design matrix for limma

group <- factor(paste(metadata$cell_line, metadata$condition, sep="."))
design <- model.matrix(~ 0 + condition, data = metadata)

# Normalization and Voom Transformation (for RNA-seq data):

## Create DGEList object and normalize
dge <- DGEList(counts = raw)
dge <- calcNormFactors(dge)

## Transform counts using voom = make raw counts to log CPM
v <- voom(dge, design, plot = TRUE)

## make cell line blocking term
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

write.csv(results_KEAP1KO, "1.1_limma_new/DE_limma_250530_blockingterm.csv")

# ## MRTX1133 vs. Control
# 
# # Define contrast matrix
# 
# contrast.matrix_MRTX1133 <- makeContrasts(
#   KO_vs_Control = (Pa14C.MRTX1133 + Pa16C.MRTX1133 + PANC1.MRTX1133 + SW1990.MRTX1133)/4 -
#     (Pa14C.control + Pa16C.control + PANC1.control + SW1990.control)/4,
#   levels=design)
# 
# 
# # Apply Contrasts and Compute Statistics:
# 
# fit_MRTX1133 <- contrasts.fit(fit, contrast.matrix_MRTX1133)
# fit_MRTX1133 <- eBayes(fit_MRTX1133)
# 
# # Extract result
# 
# results_MRTX1133  <- topTable(fit_MRTX1133, adjust="BH", number=Inf)
# 
# write.csv(results_MRTX1133, "DE_limma_MRTX1133_241027.csv")
# 
# ## KEAP1KO_MRTX1133 vs. Control
# 
# # Define contrast matrix
# 
# contrast.matrix_KEAP1KO_MRTX1133 <- makeContrasts(
#   KO.MRTX1133_vs_Control = (Pa14C.KEAP1KO_MRTX1133 + Pa16C.KEAP1KO_MRTX1133 + PANC1.KEAP1KO_MRTX1133 + SW1990.KEAP1KO_MRTX1133)/4 -
#     (Pa14C.control + Pa16C.control + PANC1.control + SW1990.control)/4,
#   levels=design)
# 
# 
# # Apply Contrasts and Compute Statistics:
# 
# fit_KEAP1KO_MRTX1133 <- contrasts.fit(fit, contrast.matrix_KEAP1KO_MRTX1133)
# fit_KEAP1KO_MRTX1133 <- eBayes(fit_KEAP1KO_MRTX1133)
# 
# # Extract result
# 
# results_KEAP1KO_MRTX1133   <- topTable(fit_KEAP1KO_MRTX1133, adjust="BH", number=Inf)
# 
# write.csv(results_KEAP1KO_MRTX1133 , "DE_limma_KEAP1KO_MRTX1133_241027.csv")
# 
# 

################## volcano plots ######################
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggrepel)


NRF2_up <- read.csv("3_NRF2up/KNR_signature_240612.csv") %>% mutate(X=NULL) 
NRF2_up.ls <- NRF2_up$ensembl_gene_id

gene_info <- read.csv("ref/gene_info.csv") %>% mutate(X = NULL)

# scaled log FC
temp <- read.csv("1.1_limma_new/DE_limma_250530_blockingterm.csv") %>%
  select(-X.1)%>%
  mutate(trend = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ "UP",
    logFC <0 & adj.P.Val < 0.05 ~ "DN",
    adj.P.Val > 0.05 ~"NOT_SIG"
  ), logp = -log10(adj.P.Val)) %>%
  mutate(scaled_p=(logp-min(logp))/(max(logp)-min(logp))) %>%
  mutate(scaled_fc=scaled_p*logFC) 


# temp <- results_KEAP1KO %>%
#   mutate(trend = case_when(
#     logFC > 0 & adj.P.Val < 0.05 ~ "UP",
#     logFC <0 & adj.P.Val < 0.05 ~ "DN",
#     adj.P.Val > 0.05 ~"NOT_SIG"
#     )) %>%
#   mutate(scaled_p=(logp-min(logp))/(max(logp)-min(logp))) %>%
#   mutate(scaled_fc=scaled_p*logFC)
  

temp_count <- temp %>% count(trend)


KNR <- temp %>%
  filter(ensembl_gene_id %in% NRF2_up$ensembl_gene_id) 

KEAP1 <-temp %>%
  filter(ensembl_gene_id == "ENSG00000079999")


p1 <- ggplot(temp) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
  geom_point(data = KNR, aes(x = logFC, y = -log10(adj.P.Val)), fill = "red", color = "black",alpha = 0.9, size = 3, shape = 21) +
  geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
  # geom_label_repel(data=KEAP1,
  #                  aes(x = logFC, y = -log10(adj.P.Val), label = "KEAP1"),
  #                  max.overlaps = Inf, show.legend = FALSE,point.padding = 0.5,
  #                  box.padding=1.2, size=5, color = "white", fill = "purple3",
  #                  segment.color = "purple", segment.size= 1.5)+
  xlim(-5,5)+
  ylim(-0.5,20)+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0) +
  labs(y="-log(adj. p value)", x="log2FC")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y   = element_text(size = 15),
        axis.text.x   = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        axis.title.x  = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.position = "none"
  )

print(p1)

#Save photos

png(paste0("1.1_limma_new/volcano_newKEAP1KO_v_control_overlap_PKD.png"), width = 2000, height = 1600, res = 300)

print(p1)

dev.off()

### determine a new PDAC KEAP1-deficiency signature
### top 200 scaled_fc, p values needs to be smaller than 0.05

sig.data <- temp %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(Rank = rank(desc(scaled_fc)))

sig.data_UP <- sig.data %>%
  filter(Rank <=200) %>%
  mutate(Signature = "PKD_UP.new")

sig.data_UP.ls <- sig.data_UP$ensembl_gene_id

sig.data_DN <- sig.data %>%
  filter(Rank > (max(Rank)-200))%>%
  mutate(Signature = "PKD_DN.new")

sig.data_DN.ls <- sig.data_DN$ensembl_gene_id

new_old_overlap.ls <- intersect(NRF2_up.ls, sig.data_UP.ls) # 67 genes overlapp with the old signature



# singscore
library(singscore)

G12Ci <- read.csv("7_Singscore/G12Ci_patient_data_240329.csv") %>%
  column_to_rownames("X")

ranked_data <- rankGenes(G12Ci)

KNRscore_G12Ci <- simpleScore(ranked_data,
                              upSet = sig.data_UP.ls)


KNR_score_result_wide <- KNRscore_G12Ci %>%
  select(sample_ID, TotalScore)mutate(suffix = sub(".*\\.", "", sample_ID),
         sample_ID = sub("\\.[^.]*$", "", sample_ID))

write.csv(KNRscore_G12Ci, "1.1_limma_new/G12Ci_patient_limma_new_signature2.csv" )
# 
# 
# 
# 
# #################### 241115, interaction term_KEAP1KO #####################
# 
# # Load necessary libraries
# library(limma)
# library(edgeR)
# library(tidyverse)
# 
# # Load data
# raw <- read.csv("KEAP1KO_counts_noNA_2.csv") %>%
#   mutate(X = NULL) %>%
#   column_to_rownames("ensembl_gene_id") 
# 
# metadata <- read.csv("1_limma/metadata_limma.csv")
# 
# # Ensure that the samples in counts data and metadata are aligned
# if (!all(colnames(raw) == metadata$sample)) {
#   # Reorder the columns of raw to match the order of metadata$sample
#   raw <- raw[, metadata$sample]
# }
# 
# # Create factors for cell line and condition
# cell_line <- factor(metadata$cell_line)
# condition <- factor(metadata$condition)
# 
# # Create the design matrix including cell line and condition as factors
# design <- model.matrix(~ cell_line + condition)
# 
# # Normalization and voom transformation
# dge <- DGEList(counts = raw)
# dge <- calcNormFactors(dge)
# v <- voom(dge, design, plot = TRUE)
# 
# # Fit the linear model
# fit <- lmFit(v, design)
# fit <- eBayes(fit)
# 
# # Identify the coefficient corresponding to the treatment effect
# colnames(fit$coefficients)  # Check the coefficient names
# coef_name <- "conditionKEAP1KO"  # Update based on your coefficient names
# 
# # Extract differential expression results for the treatment effect
# results <- topTable(fit, coef = coef_name, adjust = "BH", number = Inf)
# 
# write.csv(results, "DE_limma_interacting_term_241115.csv")
# 
# # View top differentially expressed genes
# head(results)
# 
# #volcano plots
# 
# NRF2_up <- read.csv("3_NRF2up/KNR_signature_240612.csv") %>% mutate(X=NULL) 
# gene_info <- read.csv("ref/gene_info.csv") %>% mutate(X = NULL)
# 
# temp <- read.csv("1_limma/DE_limma_interacting_term_241115.csv") %>%
#   rename(ensembl_gene_id = X)
# 
# 
# KNR <- temp %>%
#   filter(ensembl_gene_id %in% NRF2_up$ensembl_gene_id) 
# 
# KEAP1 <-temp %>%
#   filter(ensembl_gene_id == "ENSG00000079999")
# 
# 
# p1 <- ggplot(temp) +
#   geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
#   geom_point(data = KNR, aes(x = logFC, y = -log10(adj.P.Val)), fill = "red", color = "black",alpha = 0.9, size = 3, shape = 21) +
#   geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
#   # geom_label_repel(data=KEAP1,
#   #                  aes(x = logFC, y = -log10(adj.P.Val), label = "KEAP1"),
#   #                  max.overlaps = Inf, show.legend = FALSE,point.padding = 0.5,
#   #                  box.padding=1.2, size=5, color = "white", fill = "purple3",
#   #                  segment.color = "purple", segment.size= 1.5)+
#   xlim(-5,5)+
#   ylim(-0.5,40)+
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   geom_vline(xintercept = 0) +
#   labs(y="-log(adj. p value)", x="log2FC")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y   = element_text(size = 15),
#         axis.text.x   = element_text(size = 15),
#         axis.title.y  = element_text(size = 15),
#         axis.title.x  = element_text(size = 15),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         legend.position = "none"
#   )
# 
# print(p1)
# 
# #Save photos
# png(paste0("1_limma/limma_consolidated_volcano_factor_cells_KEAP1KO_241116.png"), width = 2000, height = 1600, res = 300)
# 
# print(p1)
# 
# dev.off()
# 
# 
# 
# 
# ################ KEAP1KO+ MRTX1133, 241217############################
# # Clear the workspace
# rm(list = ls())
# 
# # Load necessary libraries
# library(limma)
# library(edgeR)
# library(tidyverse)
# library(ggplot2)
# 
# # -----------------------------
# # 1. Load and Prepare the Data
# # -----------------------------
# 
# # Load count data
# raw <- read.csv("All_counts_noNA_2.csv") %>%
#   select(-X) %>%  # Remove the 'X' column if present
#   column_to_rownames("ensembl_gene_id") 
# 
# # Load metadata
# metadata <- read.csv("1_limma/metadata_limma_all.csv")
# 
# # Ensure that the samples in counts data and metadata are aligned
# if (!all(colnames(raw) == metadata$sample)) {
#   # Reorder the columns of raw to match the order of metadata$sample
#   raw <- raw[, metadata$sample]
# }
# 
# 
# 
# # Create factors for cell line and conditions
# cell_line <- factor(metadata$cell_line)
# knockout <- factor(metadata$knockout)
# treatment <- factor(metadata$treatment)
# 
# 
# # Combine knockout and treatment into a single group factor
# metadata <- metadata %>%
#   mutate(group = factor(paste(knockout, treatment, sep = "_")))
# 
# group <- factor(metadata$group)
# 
# # Verify the group levels
# print("Group Levels:")
# print(levels(metadata$group))
# 
# 
# # Create the design matrix including cell line and group
# design <- model.matrix(~ 0 + cell_line + group)
# 
# # Optional: View the design matrix column names to verify
# print("Design Matrix Columns:")
# print(colnames(design))
# 
# # Create a DGEList object
# dge <- DGEList(counts = raw)
# 
# # Calculate normalization factors
# dge <- calcNormFactors(dge)
# 
# # Apply voom transformation
# v <- voom(dge, design, plot = TRUE)
# 
# # Fit the linear model
# fit <- lmFit(v, design)
# 
# # Apply empirical Bayes moderation
# fit <- eBayes(fit)
# 
# 
# # Identify the coefficient corresponding to the treatment effect
# colnames(fit$coefficients)  # Check the coefficient names
# coef_name <- "groupKEAP1KO_MRTX1133"  # Update based on your coefficient names
# 
# # Extract differential expression results for the treatment effect
# results <- topTable(fit, coef = coef_name, adjust = "BH", number = Inf)
# 
# # Save the results to a CSV file
# write.csv(results, "1_limma/DE_limma_KEAP1KO_MRTX1133_vs_DMSO_control_limma_241217.csv", row.names = TRUE)
# 
# # View top differentially expressed genes
# head(results)
# 
# #volcano plots
# 
# NRF2_up <- read.csv("3_NRF2up/KNR_signature_240612.csv") %>% mutate(X=NULL) 
# gene_info <- read.csv("ref/gene_info.csv") %>% mutate(X = NULL)
# 
# # temp <- read.csv("1_limma/DE_limma_interacting_term_241115.csv") %>%
# #   rename(ensembl_gene_id = X)
# 
# temp <- results %>%
#   mutate(ensembl_gene_id = rownames(results))
# 
# 
# KNR <- temp %>%
#   filter(ensembl_gene_id %in% NRF2_up$ensembl_gene_id) 
# 
# KEAP1 <-temp %>%
#   filter(ensembl_gene_id == "ENSG00000079999")
# 
# 
# p1 <- ggplot(temp) +
#   geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
#   geom_point(data = KNR, aes(x = logFC, y = -log10(adj.P.Val)), fill = "red", color = "black",alpha = 0.9, size = 3, shape = 21) +
#   geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
#   # geom_label_repel(data=KEAP1,
#   #                  aes(x = logFC, y = -log10(adj.P.Val), label = "KEAP1"),
#   #                  max.overlaps = Inf, show.legend = FALSE,point.padding = 0.5,
#   #                  box.padding=1.2, size=5, color = "white", fill = "purple3",
#   #                  segment.color = "purple", segment.size= 1.5)+
#   xlim(-5,5)+
#   ylim(-0.5,40)+
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   geom_vline(xintercept = 0) +
#   labs(y="-log(adj. p value)", x="log2FC")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y   = element_text(size = 15),
#         axis.text.x   = element_text(size = 15),
#         axis.title.y  = element_text(size = 15),
#         axis.title.x  = element_text(size = 15),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         legend.position = "none"
#   )
# 
# print(p1)
# 
# #Save photos
# png(paste0("1_limma/limma_consolidated_volcano_factor_cells_MRTX1133+KEAP1_KO_241217.png"), width = 2000, height = 1600, res = 300)
# 
# print(p1)
# 
# dev.off()
# 
# 
# ### verification
# 
# KRAS_up <- read.csv("Ref/KRAS_up.csv") %>%
#   mutate(ensembl_gene_id = gsub("\\..*", "", ensembl_gene_id)) %>%
#   left_join(temp, by = "ensembl_gene_id")
# 
# KRAS_dn <- read.csv("Ref/KRAS_dn.csv") %>%
#   mutate(ensembl_gene_id = gsub("\\..*", "", ensembl_gene_id)) %>%
#   left_join(temp, by = "ensembl_gene_id")
# 
# test <- ggplot(temp) +
#   geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
#   geom_point(data = KRAS_up, aes(x = logFC, y = -log10(adj.P.Val)), fill = "green", color = "black",alpha = 0.9, size = 3, shape = 21) +
#   geom_point(data = KRAS_dn, aes(x = logFC, y = -log10(adj.P.Val)), fill = "tomato3", color = "black",alpha = 0.9, size = 3, shape = 21) +
#   geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
#   xlim(-5,5)+
#   ylim(-0.5,40)+
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   geom_vline(xintercept = 0) +
#   labs(y="-log(adj. p value)", x="log2FC", fill = "Gene Groups")+
#   # Apply classic theme and adjust legend position
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.text.y = element_text(size = 15),
#     axis.text.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.title.x = element_text(size = 15),
#     panel.border = element_rect(colour = "black", fill = NA, size = 2),
#     legend.position = "right"  # Position of the legend
#   )
# 
# print(test)
# 
# #Save photos
# png(paste0("1_limma/limma_consolidated_volcano_factor_cells_MRTX1133+KEAP1_KO_KRAS_241217.png"), width = 2000, height = 1600, res = 300)
# 
# print(p1)
# 
# dev.off()
# 
# 
# #### plotting SINGH genes
# 
# ### KEAP1KO vs. control
# gene_info <- read.csv("ref/gene_info.csv") %>% mutate(X = NULL)
# 
# temp <- read.csv("1_limma/DE_limma_interacting_term_241115.csv") %>%
#   rename(ensembl_gene_id = X)
# SINGH <- read.table("3_NRF2up/signatures/SINGH_2021.txt") %>% 
#   left_join(gene_info, by = c("V1" = "Gene")) %>%
#   left_join(temp, by = "ensembl_gene_id")
# 
# KEAP1 <-temp %>%
#   filter(ensembl_gene_id == "ENSG00000079999")
# 
# 
# p2 <- ggplot(temp) +
#   geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
#   geom_point(data = SINGH, aes(x = logFC, y = -log10(adj.P.Val)), fill = "green2", color = "black",alpha = 0.9, size = 3, shape = 21) +
#   geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
#   # geom_label_repel(data=KEAP1,
#   #                  aes(x = logFC, y = -log10(adj.P.Val), label = "KEAP1"),
#   #                  max.overlaps = Inf, show.legend = FALSE,point.padding = 0.5,
#   #                  box.padding=1.2, size=5, color = "white", fill = "purple3",
#   #                  segment.color = "purple", segment.size= 1.5)+
#   xlim(-5,5)+
#   ylim(-0.5,40)+
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   geom_vline(xintercept = 0) +
#   labs(y="-log(adj. p value)", x="log2FC")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y   = element_text(size = 15),
#         axis.text.x   = element_text(size = 15),
#         axis.title.y  = element_text(size = 15),
#         axis.title.x  = element_text(size = 15),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         legend.position = "none"
#   )
# 
# print(p2)
# 
# png(paste0("1_limma/limma_consolidated_volcano_factor_cells_KEAP1KO_Singh_241217.png"), width = 2000, height = 1600, res = 300)
# 
# print(p2)
# 
# dev.off()
# 
# # KEAP1KO + MRTX1133
# 
# temp <- read.csv("1_limma/DE_limma_KEAP1KO_MRTX1133_vs_DMSO_control_limma_241217.csv") %>%
#   rename(ensembl_gene_id = X)
# SINGH <- read.table("3_NRF2up/signatures/SINGH_2021.txt") %>% 
#   left_join(gene_info, by = c("V1" = "Gene")) %>%
#   left_join(temp, by = "ensembl_gene_id")
# 
# KEAP1 <-temp %>%
#   filter(ensembl_gene_id == "ENSG00000079999")
# 
# 
# p3 <- ggplot(temp) +
#   geom_point(aes(x = logFC, y = -log10(adj.P.Val)), fill = "grey", color = "black", alpha = 0.5, size = 3, shape = 21)+
#   geom_point(data = SINGH, aes(x = logFC, y = -log10(adj.P.Val)), fill = "green2", color = "black",alpha = 0.9, size = 3, shape = 21) +
#   geom_point(data = KEAP1, aes(x = logFC, y = -log10(adj.P.Val)), color = "purple3", size = 3) +
#   # geom_label_repel(data=KEAP1,
#   #                  aes(x = logFC, y = -log10(adj.P.Val), label = "KEAP1"),
#   #                  max.overlaps = Inf, show.legend = FALSE,point.padding = 0.5,
#   #                  box.padding=1.2, size=5, color = "white", fill = "purple3",
#   #                  segment.color = "purple", segment.size= 1.5)+
#   xlim(-5,5)+
#   ylim(-0.5,40)+
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   geom_vline(xintercept = 0) +
#   labs(y="-log(adj. p value)", x="log2FC")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y   = element_text(size = 15),
#         axis.text.x   = element_text(size = 15),
#         axis.title.y  = element_text(size = 15),
#         axis.title.x  = element_text(size = 15),
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         legend.position = "none"
#   )
# 
# print(p3)
# 
# png(paste0("1_limma/limma_consolidated_volcano_factor_cells_KEAP1KO+MRTX1133_Singh_241217.png"), width = 2000, height = 1600, res = 300)
# 
# print(p2)
# 
# dev.off()
# 
# 
# ### look at some numbers
# 
# KEAP1_KO <- read.csv("1_limma/DE_limma_interacting_term_241115.csv")
# 
# KEAP1_KO_up <- KEAP1_KO %>%
#   filter(P.Value >0.05 & logFC >0)
# 
# KEAP1_KO_dn <- KEAP1_KO %>%
#   filter(P.Value >0.05 & logFC <0)
