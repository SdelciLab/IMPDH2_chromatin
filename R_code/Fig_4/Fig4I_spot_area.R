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
    Row %in% 3 ~ "KO_empty"
  )
  )
