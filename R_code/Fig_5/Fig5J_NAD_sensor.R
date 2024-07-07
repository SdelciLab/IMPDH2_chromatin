library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(pracma)


# NAD Ratio ---------------------------------------------------------------
# 2h raw file
data_recon_2h_raw <- read.delim("/Users/alisc/Downloads/IMPDH2/IMPDH2/Data/Objects_Population - Nuclei Selected_reconssensor_2h.txt", skip = 9)


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
 # rename("NuclearArea" = "Nuclei.Selected...nuclear_areainum2n.Area..µm..",
  #       "Intensity488" = "Nuclei.Selected...Intensity.Nucleus.Alexa.488.Mean",
   #      "Intensity555" = "Nuclei.Selected...Intensity.Nucleus.Alexa.555.Mean",
  #       "YFP_mCherry" = "Nuclei.Selected...yfp.mcherry") %>% 

  filter(Nuclei.Selected...nuclear_areainum2n.Area..µm.. < 250) %>%
  # count number of samples
  group_by(Cell, Treatment) %>%
  add_tally() 

# DF containing data after removing outliers
data_recon_2h_remove_outlier <- data_recon_2h %>%
  
  group_by(Cell, Treatment) %>%
  mutate(median_yfp_mcherry = median(Nuclei.Selected...yfp.mcherry),
         sd_yfp_mcherry = sd(Nuclei.Selected...yfp.mcherry)) %>%
  
  filter(Nuclei.Selected...yfp.mcherry >= median_yfp_mcherry - (3 * sd_yfp_mcherry) &
           Nuclei.Selected...yfp.mcherry <= median_yfp_mcherry + (3 * sd_yfp_mcherry)) %>%
  ungroup() 


# Print plots to PDF
pdf("/Users/alisc/Desktop/CRG/Operetta/NAD_sensor/NAD_sensor_stats.pdf")
ggplot(data_recon_2h_remove_outlier %>% 
         filter(Cell != "WT_non-transfected"), 
       aes(Cell, Nuclei.Selected...yfp.mcherry, fill = Treatment)) + 
 # geom_point(aes(colour = Treatment),
  #           alpha=0.1, 
   #          position = position_jitterdodge(dodge.width=0.8)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette="Accent") +
  scale_colour_brewer(palette="Accent") +
  stat_compare_means(method = "wilcox", aes(group = Treatment), label = "p.signif") +
  labs(title="2h; 3SD removed;\nShowing all data") +
  theme_classic()+
  ylab("ATP/NAD ratio")+
  xlab("")
dev.off()
data_recon_2h_remove_outlier$Identifier <- paste(data_recon_2h_remove_outlier$Cell, data_recon_2h_remove_outlier$Treatment)
stats <- compare_means(data = data_recon_2h_remove_outlier, Nuclei.Selected...yfp.mcherry ~ Identifier, method = "wilcox.test")
#%>% 
 # slice(28,32)

#cell cycle analysis
data_recon_2h <- read.delim("/Users/alisc/Desktop/CRG/Operetta/NAD_sensor/Objects_Population - Nuclei Selected Selected_2h_hoechst.txt", skip = 9)

data_recon_2h <- data_recon_2h %>%
  mutate(Cell = case_when(
    Row == 2 ~ "WT_empty",
    Row == 3 ~ "KO_empty",
    Row == 4 ~ "KO-WT",
    Row == 5 ~ "KO-NLS",
    Row == 6 ~ "WT_non-transfected"
  )
  )
data_recon_2h <- data_recon_2h[data_recon_2h$Cell != "WT_non-transfected",]
data_recon_2h <- data_recon_2h %>%
  mutate(Treatment = case_when(
    Column %in% c(4,5,6) ~ "DMSO",
    Column %in% c(7,8,9) ~ "Eto_10uM"
  )
  )

#exclude cells that are too big
data_recon_2h_filtered <- data_recon_2h[data_recon_2h$Nuclei.Selected.Selected...nuclear_areainum2n.Area..µm..< 400, ]
ggplot(data_recon_2h_filtered, aes(Nuclei.Selected.Selected...nuclear_areainum2n.Area..µm..)) + geom_density()

#integrated hoechst signal
data_recon_2h_filtered$Hoechst_corrected <- data_recon_2h_filtered$Nuclei.Selected.Selected...Intensity.Nucleus.HOECHST.33342.Mean - data_recon_2h_filtered$Nuclei.Selected.Selected...Intensity.Ring.Region.HOECHST.33342.Mean

data_recon_2h_filtered$ratio <- data_recon_2h_filtered$Hoechst_corrected*data_recon_2h_filtered$Nuclei.Selected.Selected...nuclear_areainum2n.Area..µm..
data_recon_2h_filtered <- data_recon_2h_filtered[data_recon_2h_filtered$ratio<1e+07,]

data_2h_DMSO <- data_recon_2h_filtered[data_recon_2h_filtered$Treatment == "DMSO",]
data_2h_Etoposide <- data_recon_2h_filtered[data_recon_2h_filtered$Treatment == "Eto_10uM",]

ggplot(data_2h_DMSO, aes(ratio, color=factor(Column))) +geom_density() + facet_wrap(~Cell)
ggplot(data_2h_Etoposide, aes(ratio, color=factor(Column))) +geom_density() + facet_wrap(~Cell)

#normalise G1 peaks across all replicates
data_2h_DMSO$Identifier <- paste(data_2h_DMSO$Cell, data_2h_DMSO$Column)
a <- ggplot(data_2h_DMSO, aes(ratio, color=Identifier)) +geom_density()
a

density_data <- ggplot_build(a)$data[[1]]
group_levels <- levels(factor(data_2h_DMSO$Identifier))
density_data$Identifier <- group_levels[density_data$group]

#check for peak value per condtion
density <- density_data %>%
  group_by(Identifier) %>%
  mutate(MaxY = max(y)) %>% 
  filter(y == MaxY) %>% 
  summarise(x = first(x), MaxY = first(MaxY), .groups = 'drop')

#check for G1 peak and normalise
# Shift for WT replicates 
density_WT <- density[density$Identifier %in% c("WT_empty 4", "WT_empty 5", "WT_empty 6"),]
wt_empty_peak_x <- min(density_WT$x)
density_WT$shift = wt_empty_peak_x - density_WT$x

# Shift for KO replicates 
density_KO <- density[density$Identifier %in% c("KO_empty 4", "KO_empty 5", "KO_empty 6"),]
ko_empty_peak_x <- min(density_KO$x)
density_KO$shift = ko_empty_peak_x - density_KO$x

# Shift for KO-NLS replicates 
density_KONLS <- density[density$Identifier %in% c("KO-NLS 4", "KO-NLS 5", "KO-NLS 6"),]
KONLS_peak_x <- min(density_KONLS$x)
density_KONLS$shift = KONLS_peak_x - density_KONLS$x

# Shift for KO-WT replicates 
density_KOWT <- density[density$Identifier %in% c("KO-WT 4", "KO-WT 5", "KO-WT 6"),]
kowt_peak_x <- min(density_KOWT$x)
density_KOWT$shift = kowt_peak_x - density_KOWT$x

shift_data <- rbind(density_KOWT, density_KONLS, density_KO, density_WT)

# Join data with shifts
normalised_data_repicates <- data_2h_DMSO %>%
  left_join(shift_data, by = "Identifier") %>%
  mutate(adjusted_value = ratio + shift) 

#plot normalised 
ggplot(normalised_data_repicates, aes(adjusted_value, color=factor(Column))) +geom_density() + facet_wrap(~Cell)
#plot all replicates together
b <- ggplot(normalised_data_repicates, aes(adjusted_value, color=Cell)) +geom_density()
b
#normalise Wt empty G1 peak
density_data <- ggplot_build(b)$data[[1]]
group_levels <- levels(factor(normalised_data_repicates$Cell))
density_data$Cell <- group_levels[density_data$group]

#check for peak value per condition
density <- density_data %>%
  group_by(Cell) %>%
  mutate(MaxY = max(y)) %>% 
  filter(y == MaxY) %>% 
  summarise(x = first(x), MaxY = first(MaxY), .groups = 'drop')

# WT peak x-coordinate
WT_empty_peak_x <- density$x[density$Cell == "WT_empty"]

# Calculate shifts needed for each condition
density$shift2 = WT_empty_peak_x - density$x

# Join data with shifts
normalised_data_G1 <- normalised_data_repicates %>%
  left_join(density, by = "Cell") %>%
  mutate(adjusted_value_G1 = adjusted_value + shift2) 

#define gates
gate1 <- WT_empty_peak_x -400000
gate2 <- WT_empty_peak_x[1] + 400000

#find G2 peak
data_WT_empty <- normalised_data_G1 %>%
  dplyr::filter(Cell %in% c("WT_empty"))
dens <- density(data_WT_empty$adjusted_value_G1)
peaks <- findpeaks(dens$y)
peaks <- as.data.frame(peaks)
peak_x_values <- dens$x[peaks$V2]

G2peak <- 4022573
gate3 <- G2peak -450000
gate4 <- G2peak + 450000

ggplot(normalised_data_G1, aes(adjusted_value_G1, color=Cell)) + geom_density() +  # Density curve
  geom_vline(xintercept = WT_empty_peak_x, col = "red", linetype = "dashed") +
  geom_vline(xintercept = G2peak, col = "red", linetype = "dashed") +
  geom_vline(xintercept = gate1, col = "blue") +
  geom_vline(xintercept = gate2, col = "blue") + geom_vline(xintercept = gate3, col = "blue") +
  geom_vline(xintercept = gate4, col = "blue") + theme_bw()

normalised_data_subset_2h_DMSO <- normalised_data_G1 %>%
  mutate(Cell_cycle = case_when(
    adjusted_value_G1 >= gate1 & adjusted_value_G1 <= gate2 ~ "G1",
    adjusted_value_G1 > gate2 & adjusted_value_G1 < gate3 ~ "S",
    adjusted_value_G1 >= gate3 & adjusted_value_G1 <= gate4 ~ "G2/M",
    TRUE ~ "X"
  ))
normalised_data_subset_2h_DMSO <- normalised_data_subset_2h_DMSO[normalised_data_subset_2h_DMSO$Cell_cycle != "X",]

my_comp2 <-  list(c("WT_empty", "KO_empty"),c("KO-WT","KO-NLS"))
ggplot(normalised_data_subset_2h_DMSO, aes(Cell, Nuclei.Selected.Selected...yfp.mcherry, fill= Cell_cycle)) +  geom_boxplot() + scale_color_brewer(palette="Accent")+ labs(title="DMSO") +theme_classic() 

data_recon_2h_filtered_signal_outlier_removed <- normalised_data_subset_2h_DMSO %>%
  group_by(Cell, Cell_cycle) %>%
  mutate(median_yfp_mcherry = median(Nuclei.Selected.Selected...yfp.mcherry),
         sd_yfp_mcherry = sd(Nuclei.Selected.Selected...yfp.mcherry)) %>%
  ungroup() %>%
  filter(Nuclei.Selected.Selected...yfp.mcherry >= median_yfp_mcherry - 3 * sd_yfp_mcherry &
           Nuclei.Selected.Selected...yfp.mcherry <= median_yfp_mcherry + 3 * sd_yfp_mcherry)

ggplot(data_recon_2h_filtered_signal_outlier_removed, aes(Cell, Nuclei.Selected.Selected...yfp.mcherry, fill= Cell_cycle)) +  geom_boxplot() + scale_fill_brewer(palette="Accent")+ labs(title="DMSO") +theme_classic() 

summary_data <- data_recon_2h_filtered_signal_outlier_removed %>%
  group_by(Cell, Cell_cycle) %>%
  summarize(Count = n(), .groups = 'drop') 
