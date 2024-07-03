library(tidyverse)
library(ggpubr)
library(ggplot2)
#Figure 4H column of interest: Cyto Selected - Intensity Cytoplasm Alexa 488 Mean
#WT: B6, C6, D6; KO: E6, F6, G6, file Objects_Population - Cyto Selected_Fig4H

#load the data
data <- read.delim("/Users/alisc/Downloads/Objects_Population - Cyto Selected_Fig4H.txt", skip = 9)

data <- data[data$Column == 6,]
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
  mutate(median_1 = median(Cyto.Selected...Intensity.Cytoplasm.Alexa.488.Mean, na.rm= T),
         sd_1 = sd(Cyto.Selected...Intensity.Cytoplasm.Alexa.488.Mean, na.rm=T)) %>%
  
  filter(Cyto.Selected...Intensity.Cytoplasm.Alexa.488.Mean >= median_1 - (3 * sd_1) &
           Cyto.Selected...Intensity.Cytoplasm.Alexa.488.Mean <= median_1 + (3 * sd_1)) %>%
  ungroup() 
ordered <- c("WT", "KO")
data_no_outlier$Cell <- factor(data_no_outlier$Cell, levels = ordered)

pdf("/Users/alisc/Desktop/CRG/final_figures/Fig4H_PARP1_spot_intensity.pdf")
ggplot(data_no_outlier, aes(x = Cell, y = Cyto.Selected...Intensity.Cytoplasm.Alexa.488.Mean, fill=Cell)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "PARP1 cytosolic intensity") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme(legend.position = "none")  + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  scale_fill_brewer(palette="Accent")
dev.off()

table(data_no_outlier$Cell)
