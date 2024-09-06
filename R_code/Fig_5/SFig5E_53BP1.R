library(tidyverse)
library(ggpubr)
library(ggplot2)

#Row/Column of interest
#E2F2G2 WT
#E3F3G3 KO
#E4F4G4 KO-WT
#E6F6G6 KO-NLS
data <- read.delim("/Users/alisc/Desktop/CRG/Operetta/MDAMB231 WT, KO, KO-WT, KO-NLS/Objects_Population - full nuclei.txt", sep="\t",  skip = 9)
row_noG <- c(5,6,7)
data_noG <- data[data$Row %in% row_noG,]
data_cond <- data_noG %>%
  mutate(Cell = case_when(
    Column == 2 ~ "WT",
    Column == 3 ~ "KO",
    Column == 4 ~ "KO-WT",
    Column == 6 ~ "KO-NLS"))
data_cond <- data_cond %>% drop_na(Cell)
data_cond$Cell <- factor(data_cond$Cell, levels = c("WT", "KO","KO-WT","KO-NLS"))

# comparisons
my_comp <- list(c("WT", "KO"), c("KO", "KO-WT"), c("KO-WT", "KO-NLS"),c("KO", "KO-NLS"))

ggplot(data_cond, aes(Cell, full.nuclei...Number.of.Spots..2. , fill=Cell)) +
  geom_boxplot(width = 0.75, fatten = 2) + 
  theme_minimal() +  
  labs(y = "53BP1 spots") + 
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) 

#remove outlier
data_cond_1 <- data_cond %>%
  group_by(Cell) %>%
  mutate(median_53BP1 = median(full.nuclei...Number.of.Spots..2.),
         sd_53BP1 = sd(full.nuclei...Number.of.Spots..2.)) %>%
  filter(full.nuclei...Number.of.Spots..2. >= median_53BP1 - (3 * sd_53BP1) &
           full.nuclei...Number.of.Spots..2. <= median_53BP1 + (3 * sd_53BP1)) %>%
  ungroup()

ggplot(data_cond_1, aes(Cell, full.nuclei...Number.of.Spots..2. , fill=Cell)) +
  geom_boxplot(width = 0.75, fatten = 2) + 
  theme_minimal() +  
  labs(y = "53BP1 spots") + 
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01)  +
  scale_fill_brewer(palette="Accent")
