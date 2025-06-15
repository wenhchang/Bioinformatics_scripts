#### 2.5 visualize GSEA analysis for DE results of each cell line

rm(list = ls()) #clear the environment

# load necessary packages
library(tidyverse)

# load data
GSEA_H <- read.csv("GSEA_H.csv")

## clean up data 
GSEA_H_data <- GSEA_H_data %>%
  select(pathway, padj, NES, source) %>%
  mutate(source = gsub("Pa14C_new", "Pa14C", source), 
         source = gsub("PANC1", "PANC-1", source), 
         pathway = gsub("HALLMARK_", "", pathway))

# set the order 
## manually set the order of hallmark pathway in the figure
## for the order, I rank the gene sets from the highest medium NES value to the lowest

GSEA_H_data_sumpNES <- GSEA_H_data %>%
  select(pathway, NES, source) %>%
  pivot_wider(names_from = source, values_from = NES) %>%
  # mutate(mean = rowMed(select(., SW1990:Pa14C))) %>%
  mutate(median = apply(select(., SW1990:Pa14C), 1, median)) %>%
  arrange(desc(median)) %>%
  mutate(rank_NES = row_number())

order <- factor(GSEA_H_data_sumpNES$pathway) 

# visualization of GSEA results for Hallmark gene sets
p1 <- GSEA_H_data %>% #your data set
  ggplot(aes(x = pathway, y = source, fill = NES, size = -log10(padj))) + #I combined all four cells in the data sets by adding a "source" column
  geom_point(shape = 21) +
  theme_classic() +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "red2", midpoint = 0, 
                       name = "NES", limits = c(-2, 2), oob = scales::squish)+ #this is to manually set the color representation
  scale_x_discrete(limits = order)+  # Set x-axis order based on 'order' vector
  scale_size_continuous(range = c(0.1, 6)) + #the dot size scale
  # labs(y = "Cell Line",
  #      x = "Gene set",
  #      title = "Hallmark gene sets GSEA") +
  theme(axis.text.y = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, color = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # plot.title = element_text(size = 22), 
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.position = "bottom",
        plot.background = element_rect(fill = "transparent", color = NA), # Transparent plot background
        panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel background
        legend.background = element_rect(fill = "transparent", color = NA) # Optional: Transparent legend background
  )

#saving the plot
png("GSEA.png", width = 5500, height = 2800, res = 600, bg = "transparent")
print(p1)
dev.off()
