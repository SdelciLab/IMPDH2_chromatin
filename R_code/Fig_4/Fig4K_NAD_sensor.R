library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)


# NAD Ratio ---------------------------------------------------------------
# 2h raw file
data_recon_2h_raw <- read.delim("Data/Objects_Population - Nuclei Selected_reconssensor_2h.txt", skip = 9)

# 2h processed file
data_recon_2h <- data_recon_2h_raw %>%  
  as_tibble() %>% 
  mutate(Cell = case_when(Row == 2 ~ "WT_empty",
                          Row == 3 ~ "KO_empty",
                          Row == 4 ~ "KO_WT",
                          Row == 5 ~ "KO_NLS",
                          Row == 6 ~ "WT_non-transfected"),
         Treatment = case_when(Column %in% c(4,5,6) ~ "DMSO",
                               Column %in% c(7,8,9) ~ "Eto_10uM"),
         Cell = factor(Cell, levels = c("WT_empty","KO_empty", "KO_WT","KO_NLS", "WT_non-transfected")),
         Treatment = factor(Treatment, levels = c("DMSO","Eto_10uM"))) %>% 
  rename("NuclearArea" = "Nuclei.Selected...nuclear_areainum2n.Area..?m..",
         "Intensity488" = "Nuclei.Selected...Intensity.Nucleus.Alexa.488.Mean",
         "Intensity555" = "Nuclei.Selected...Intensity.Nucleus.Alexa.555.Mean",
         "YFP_mCherry" = "Nuclei.Selected...yfp.mcherry") %>% 

  filter(NuclearArea < 250) %>%
  # count number of samples
  group_by(Cell, Treatment) %>%
  add_tally() 

# DF containing data after removing outliers
data_recon_2h_remove_outlier <- data_recon_2h %>%
  
  group_by(Cell, Treatment) %>%
  mutate(median_yfp_mcherry = median(YFP_mCherry),
         sd_yfp_mcherry = sd(YFP_mCherry)) %>%
  
  filter(YFP_mCherry >= median_yfp_mcherry - (3 * sd_yfp_mcherry) &
           YFP_mCherry <= median_yfp_mcherry + (3 * sd_yfp_mcherry)) %>%
  ungroup() 


# Print plots to PDF
ggplot(data_recon_2h_remove_outlier %>% 
         filter(Cell != "WT_non-transfected"), 
       aes(Cell, YFP_mCherry, fill = Treatment)) + 
  geom_point(aes(colour = Treatment),
             alpha=0.1, 
             position = position_jitterdodge(dodge.width=0.8)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette="Accent") +
  scale_colour_brewer(palette="Accent") +
  labs(title="2h; 3SD removed;\nShowing all data") +
  theme_classic()+
  ylab("NAD/NAD+ ratio")+
  xlab("")


