#### RPPA analysis: 2.0 Differential expression analysis

rm(list = ls())

library(tidyverse)
library(limma)

# Load pre-processed data
dat <- readRDS("results/01_preprocessed_data.rds")
expr <- dat$expr
meta <- dat$meta


# Run differential expression analysis with limma

expr <- expr[, meta$ID] ## ensure column order

## check if the order is aligned
colnames(expr) == metadata$ID


## prepare design matrix
design <- model.matrix(~0 + sample, data = meta)
colnames(design) <- levels(meta$sample)

## limma linear model fitting
fit <- lmFit(expr, design)

## create contrast matrix: include all the contrasts of interest
contrast.matrix <- makeContrasts(
  R5W_vs_WT     = R5W - WT,
  R5W_vs_Y42C   = R5W - Y42C,
  R5W_vs_L57V   = R5W - L57V,
  Y42C_vs_WT    = Y42C - WT,
  G17E_vs_WT    = G17E - WT,
  L57V_vs_WT    = L57V - WT,
  Q63L_vs_WT    = Q63L - WT,
  WT_vs_vector  = WT - vector,
  R5W_vs_vector = R5W - vector,
  Y42C_vs_vector= Y42C - vector,
  G17E_vs_vector= G17E - vector,
  L57V_vs_vector= L57V - vector,
  Q63L_vs_vector= Q63L - vector,
  levels = design)

## run contrast fit and eBayes

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

## export contrast tables

result.table <- data.frame()

for (con in colnames(contrast.matrix)) {
  tt <- topTable(fit2, coef = con, adjust.method = "BH", number = Inf, sort.by = "logFC") %>%
    mutate(contrast = con)
  
  result.table <- rbind(tt, result.table)
  
}

write.csv(result.table, "RPPA_DE_result.csv")
