library(tidyverse)
library(ggpubr)
library(ggplot2)
#Figure 4G column of interest: Cyto Selected - Total Spot Area - Mean per Well 
#WT: B11, C11, D11; KO: E11, F11, G11, file Objects_Population - Cyto Selected_Fig4G

#load the data
data <- read.delim("/Users/alisc/Downloads/Objects_Population - Cyto Selected_Fig4G.txt", skip = 9)

data <- data[data$Column == 11,]
data_cond <- data %>%
  mutate(Cell = case_when(
    Row %in% c(2,3,4) ~ "WT",
    Row %in% c(5,6,7) ~ "KO"
  )
  )
my_comparisons <- list(c("WT", "KO"))
#remove outlier
data_no_outlier <- data_cond %>%
  
  group_by(Cell) %>%
  mutate(median_1 = median(Cyto.Selected...Total.Spot.Area),
         sd_1 = sd(Cyto.Selected...Total.Spot.Area)) %>%
  
  filter(Cyto.Selected...Total.Spot.Area >= median_1 - (3 * sd_1) &
           Cyto.Selected...Total.Spot.Area <= median_1 + (3 * sd_1)) %>%
  ungroup() 
ordered <- c("WT", "KO")
data_no_outlier$Cell <- factor(data_no_outlier$Cell, levels = ordered)
pdf("/Users/alisc/Desktop/CRG/final_figures/Fig4G_spot_area.pdf")
ggplot(data_no_outlier, aes(x = Cell, y = Cyto.Selected...Total.Spot.Area, fill=Cell)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "Area of cytosolic cleaved caspase 3 foci") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme(legend.position = "none")  + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  scale_fill_brewer(palette="Accent")
dev.off()

table(data_no_outlier$Cell)
