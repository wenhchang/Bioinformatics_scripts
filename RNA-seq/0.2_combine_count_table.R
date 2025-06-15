#### 0.2 combine results directly after STAR mapping and featureCounts to a single count table

rm(list=ls()) 

# load necessary packages

library(tidyverse)

#### make sure gene_length and the data are in the same order !!!!!!!!!

# load sample names in the count matrix
sample_name <- c("SW1990", "PANC1", "Pa16C", "Pa14C_new")


# create an emmpty dataframe for later result storage
all.df <- data.frame(ensembl_gene_id = character())

for (cell in sample_name) {

  raw <- read.csv(paste0("raw_counts/", cell,".csv")) %>%
    mutate(X=NULL) %>%
    rename(ensembl_gene_id = Geneid) %>%
    mutate(ensembl_gene_id = gsub("\\..*$", "", ensembl_gene_id)) 

  # filter out raw counts < 20
  raw_filter <- raw %>%
    mutate(average_counts = rowMeans(select(., 2:13 ), na.rm = TRUE)) %>%
    filter(average_counts >= 20) %>%
    mutate(average_counts = NULL)
  
  all.df <- all.df %>%
    full_join(raw_filter, by ="ensembl_gene_id")  


}

# write.csv(all.df, "KEAP1KO_counts.csv")
