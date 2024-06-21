library(tidyverse)
library(ggplot2)
library(ggpubr)

#Conditions
#3 cell lines: MCF7, T47D, MDAMB231
#Cells in columns: MCF7 (3), T47D (4),CAL51 (5) -> do not include, MDAMB231 (6)
# get data
data_WT_cells <- read.delim("/Users/alisc/Desktop/CRG/Operetta/Fig1B/Objects_Population - full nuclei_MCF7T47DMDAMB231.txt", skip = 9)

# select only interesting columns
data_WT_cells_int <- data_WT_cells %>%
  filter(Column %in% c(3,4,6))

# add cell line
data_WT_cells_int_cond <- data_WT_cells_int %>%
  mutate(Cell = case_when(
    Column == 3 ~ "MCF7",
    Column == 4 ~ "T47D",
    Column == 6 ~ "MDAMB231")
  )

# factor
data_WT_cells_int_cond$Cell <- factor(data_WT_cells_int_cond$Cell, levels = c("MCF7", "T47D",
                                                                              "MDAMB231"))
#sampling to represent same number of datapoints 
#MDAMB231 lowest number of datapoints with 748

# Filter data for MDAMB231, keep all data points as is
data_MDAMB231  <- data_WT_cells_int_cond %>%
  filter(Cell == "MDAMB231")

# Filter data for MCF7/T27D and then sample 748 data points
set.seed(100)
data_MCF7 <- data_WT_cells_int_cond %>%
  filter(Cell == "MCF7") %>%
  sample_n(size = 748 )

data_T47D <- data_WT_cells_int_cond %>%
  filter(Cell == "T47D") %>%
  sample_n(size = 748 )

# Combine the balanced datasets
balanced_data <- bind_rows(data_MDAMB231, data_MCF7, data_T47D)

#Plot: full nuclei - Number of Spots (H2AX spots)
# comparisons
my_comp <- list(c("MCF7", "T47D"), c("T47D", "MDAMB231"), c("MCF7", "MDAMB231"))

#pdf("/Users/alisc/Desktop/CRG/final_figures/H2AX_MCF7_T47D_MDAMB231/balanced_jitter.pdf")
ggplot(balanced_data, aes(x = Cell, y = full.nuclei...Number.of.Spots, color = Cell)) + geom_jitter(alpha=0.4,  width = 0.25) +
  geom_boxplot(color= "black", fill= NA, width = 0.5) +
  labs(x = "", y = "Number of H2AX foci per nucleus)", fill = "") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  scale_color_manual(name = "", values = c("MDAMB231" = "#354f52", 
                                           "MCF7" = "#CAD2C5",
                                           "T47D" = "#CAD2C5"))+
  theme_classic() +
  theme( legend.position = "none")
#dev.off()

