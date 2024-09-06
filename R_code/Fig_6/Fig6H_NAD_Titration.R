library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)


# NAD: Zero Mean and Unit Variance Normalization ---------------------------------------------------------------
# Read in PAR data
data_par_raw <- read.delim("/Users/alisc/Downloads/IMPDH2/IMPDH2/Data/Objects_Population - Nuclei 2.txt", skip = 9)

data_par <- data_par_raw %>%
  as_tibble() %>% 
  mutate(Cell = case_when(
    Row == 2 ~ "WT_empty",
    Row == 3 ~ "KO_empty",
    Row == 4 ~ "KO_WT",
    Row == 5 ~ "KO_NLS",
    Row == 6 ~ "WT",
    Row == 7 ~ "KO"),
    Cell = factor(Cell, levels = c("WT", "KO","WT_empty","KO_empty","KO_WT","KO_NLS"))) %>% 
  rename("NuclearArea" = "Nuclei...Nucleus.Area",
         "Intensity488" = "Nuclei...Intensity.Nucleus.Alexa.488.Mean",
         "Intensity488max" = "Nuclei...Intensity.Nucleus.Alexa.488.Maximum",
         "Intensity488min" = "Nuclei...Intensity.Nucleus.Alexa.488.Minimum") %>% 
  filter(Cell %in% c("KO_WT", "KO_NLS"),
         NuclearArea < 250) %>% 
  group_by(Cell) %>% 
  mutate(meanIntensity488 = mean(Intensity488),
         sdIntensity488 = sd(Intensity488)) %>% 
  ungroup()


# Read in titration data
data_titration_raw <- read.delim("/Users/alisc/Downloads/IMPDH2/IMPDH2/Data/Objects_Population - full nuclei.txt", skip=9)

data_titration <- data_titration_raw %>%
  as_tibble() %>% 
  mutate(NAD = case_when(
    Column == 2 & Row %in% c(5,6,7) ~ "ko-wt",
    Column == 3 ~ "1000",
    Column == 4 ~ "800",
    Column == 5 ~ "400",
    Column == 6 ~ "200",
    Column == 7 ~ "100",
    Column == 8 ~ "50",
    Column == 9 ~ "25",
    Column == 10 ~ "12.5",
    Column == 11 ~ "6.25",
    Column == 2 ~ "0"),
    NAD = factor(NAD, levels = c("ko-wt","0","6.25","12.5","25","50","100","200","400","800","1000"))) %>% 
  rename("NuclearArea" = "full.nuclei...Nucleus.Area..µm..",
         "Intensity488" = "full.nuclei...Intensity.Nucleus.Alexa.488.Mean",
         "Intensity488max" = "full.nuclei...Intensity.Nucleus.Alexa.488.Maximum",
         "Intensity488min" = "full.nuclei...Intensity.Nucleus.Alexa.488.Minimum",
         "Roundness" = "full.nuclei...Nucleus.Roundness") %>% 
  filter(Roundness > 0.8,
         NuclearArea < 300) %>% 
  group_by(NAD) %>% 
  mutate(median_488 = median(Intensity488),
         sd_488 = sd(Intensity488)) %>% 
  ungroup() %>% 
  filter(Intensity488 >= (median_488 - 5 * sd_488) & Intensity488 <= (median_488 + 5 * sd_488),
         !NAD %in% c("ko-wt","400", "1000")) %>% 
  group_by(NAD) %>% 
  mutate(meanIntensity488 = mean(Intensity488),
         sdIntensity488 = sd(Intensity488)) %>% 
  ungroup()


# Normalizing values in PAR dataset:
PAR_KO_NLS = filter(data_par, Cell=="KO_NLS") %>% 
  mutate(newIntensity488 = ((Intensity488 - meanIntensity488)/sdIntensity488)) %>% 
  select(Cell, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)

par_control_mean = PAR_KO_NLS$meanIntensity488 %>% unique()
par_control_sd = PAR_KO_NLS$sdIntensity488 %>% unique()

PAR_KO_WT = filter(data_par, Cell=="KO_WT") %>% 
  mutate(newIntensity488 = ((Intensity488 - par_control_mean)/par_control_sd)) %>% 
  select(Cell, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)

PAR_normalised <- rbind(PAR_KO_WT, PAR_KO_NLS)

# Normalizing values in Titration dataset:
Titre_0 = filter(data_titration, NAD == "0") %>% 
  mutate(newIntensity488 = ((Intensity488 - meanIntensity488)/sdIntensity488)) %>% 
  select(NAD, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)

titre_control_mean = Titre_0$meanIntensity488 %>% unique()
titre_control_sd = Titre_0$sdIntensity488 %>% unique()

normalised_NAD <- Titre_0

for(i in data_titration$NAD %>% unique()){
  
  temp_data <- filter(data_titration, NAD == i) %>% 
    mutate(newIntensity488 = ((Intensity488 - titre_control_mean)/titre_control_sd)) %>% 
    select(NAD, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)
  
  normalised_NAD <- rbind(normalised_NAD, temp_data)
  
}

normalised_NAD <- rename(normalised_NAD, "Cell" = "NAD")


# Combining the two data frames:

normalised_combined_df <- rbind(PAR_normalised, normalised_NAD) %>% 
  filter(Cell %in% c("KO_WT", "KO_NLS", "0", "12.5", "50", "100", "200", "400")) %>% 
  mutate(Cell = factor(Cell, levels = c("KO_NLS", "0", "12.5", "50", "100", "200", "400", "KO_WT"))) %>% 
  group_by(Cell) %>% 
  mutate(mean_newIntensity488 = mean(newIntensity488)) %>% 
  add_tally() %>% 
  ungroup()

#stat_normalised_NAD <- compare_means(data = normalised_combined_df, newIntensity488 ~ Cell, method = "wilcox.test") %>% 
# slice(1:7, 13, 18, 22, 25, 27:28)
stat_normalised_NAD <- compare_means(data = normalised_combined_df, newIntensity488 ~ Cell, method = "wilcox.test") %>% 
  slice(1:6, 11, 15, 18, 20, 21)

p_normalised <- ggplot(normalised_combined_df, aes(x=Cell, y=newIntensity488, fill=))+
  geom_boxplot(aes(fill = mean_newIntensity488), outlier.shape = NA) +
  coord_cartesian(ylim = c(-1.2, 6.5)) +
  scale_fill_viridis_c(direction = 1)+
  stat_pvalue_manual(stat_normalised_NAD,  
                     label = "p.signif", 
                     y.position = 2.5,
                     step.increase = 0.024,
                     tip.length = 0.01)+
  theme_bw() +
  ylab("Intensity Mono/Poly ADP-ribose")+
  xlab("Conditions")+
  ggtitle("Zero Mean and Unit Variance Normalization")+
  labs(fill='Mean y-value')

# For easier viewership, add a constant - typically equal to the minimum value.
p_normalised_constant <- ggplot(normalised_combined_df, aes(x=Cell, y = newIntensity488 + abs(min(newIntensity488))))+
  geom_boxplot(aes(fill = mean_newIntensity488 + abs(min(newIntensity488))), outlier.shape = NA) +
  scale_fill_viridis_c(direction = 1)+
  coord_cartesian(ylim = c(0, 8)) +
  stat_pvalue_manual(stat_normalised_NAD,  
                     label = "p.signif", 
                     y.position = 4,
                     step.increase = 0.024,
                     tip.length = 0.01)+
  theme_bw()+
  ylab("Intensity Mono/Poly ADP-ribose")+
  xlab("Conditions")+
  ggtitle("Zero Mean and Unit Variance Normalization (+Constant)")+
  labs(fill='Mean y-value')
p_normalised_constant


#just KO_NLS, 0, 12.5, 100, KO_WT
data_titration <- data_titration_raw %>%
  as_tibble() %>% 
  mutate(NAD = case_when(
    Column == 2 & Row %in% c(5,6,7) ~ "ko-wt",
    Column == 3 ~ "1000",
    Column == 4 ~ "800",
    Column == 5 ~ "400",
    Column == 6 ~ "200",
    Column == 7 ~ "100",
    Column == 8 ~ "50",
    Column == 9 ~ "25",
    Column == 10 ~ "12.5",
    Column == 11 ~ "6.25",
    Column == 2 ~ "0"),
    NAD = factor(NAD, levels = c("ko-wt","0","6.25","12.5","25","50","100","200","400","800","1000"))) %>% 
  rename("NuclearArea" = "full.nuclei...Nucleus.Area..µm..",
         "Intensity488" = "full.nuclei...Intensity.Nucleus.Alexa.488.Mean",
         "Intensity488max" = "full.nuclei...Intensity.Nucleus.Alexa.488.Maximum",
         "Intensity488min" = "full.nuclei...Intensity.Nucleus.Alexa.488.Minimum",
         "Roundness" = "full.nuclei...Nucleus.Roundness") %>% 
  filter(Roundness > 0.8,
         NuclearArea < 300) %>% 
  group_by(NAD) %>% 
  mutate(median_488 = median(Intensity488),
         sd_488 = sd(Intensity488)) %>% 
  ungroup() %>% 
  filter(Intensity488 >= (median_488 - 5 * sd_488) & Intensity488 <= (median_488 + 5 * sd_488),
         NAD %in% c("KO_NLS","0", "12.5", "100", "KO_WT")) %>% 
  group_by(NAD) %>% 
  mutate(meanIntensity488 = mean(Intensity488),
         sdIntensity488 = sd(Intensity488)) %>% 
  ungroup()


# Normalizing values in PAR dataset:
PAR_KO_NLS = filter(data_par, Cell=="KO_NLS") %>% 
  mutate(newIntensity488 = ((Intensity488 - meanIntensity488)/sdIntensity488)) %>% 
  select(Cell, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)

par_control_mean = PAR_KO_NLS$meanIntensity488 %>% unique()
par_control_sd = PAR_KO_NLS$sdIntensity488 %>% unique()

PAR_KO_WT = filter(data_par, Cell=="KO_WT") %>% 
  mutate(newIntensity488 = ((Intensity488 - par_control_mean)/par_control_sd)) %>% 
  select(Cell, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)

PAR_normalised <- rbind(PAR_KO_WT, PAR_KO_NLS)

# Normalizing values in Titration dataset:
Titre_0 = filter(data_titration, NAD == "0") %>% 
  mutate(newIntensity488 = ((Intensity488 - meanIntensity488)/sdIntensity488)) %>% 
  select(NAD, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)

titre_control_mean = Titre_0$meanIntensity488 %>% unique()
titre_control_sd = Titre_0$sdIntensity488 %>% unique()

normalised_NAD <- Titre_0

for(i in data_titration$NAD %>% unique()){
  
  temp_data <- filter(data_titration, NAD == i) %>% 
    mutate(newIntensity488 = ((Intensity488 - titre_control_mean)/titre_control_sd)) %>% 
    select(NAD, newIntensity488, Intensity488, meanIntensity488, sdIntensity488)
  
  normalised_NAD <- rbind(normalised_NAD, temp_data)
  
}

normalised_NAD <- rename(normalised_NAD, "Cell" = "NAD")


# Combining the two data frames:

normalised_combined_df <- rbind(PAR_normalised, normalised_NAD) %>% 
  filter(Cell %in% c("KO_WT", "KO_NLS", "0", "12.5", "50", "100", "200", "400")) %>% 
  mutate(Cell = factor(Cell, levels = c("KO_NLS", "0", "12.5", "50", "100", "200", "400", "KO_WT"))) %>% 
  group_by(Cell) %>% 
  mutate(mean_newIntensity488 = mean(newIntensity488)) %>% 
  add_tally() %>% 
  ungroup()

stat_normalised_NAD <- compare_means(data = normalised_combined_df, newIntensity488 ~ Cell, method = "wilcox.test") %>% 
  slice(1:4, 7, 9,10)

p_normalised <- ggplot(normalised_combined_df, aes(x=Cell, y=newIntensity488, fill=))+
  geom_boxplot(aes(fill = mean_newIntensity488), outlier.shape = NA) +
  coord_cartesian(ylim = c(-1.2, 6.5)) +
  scale_fill_viridis_c(direction = 1)+
  stat_pvalue_manual(stat_normalised_NAD,  
                     label = "p.signif", 
                     y.position = 3.5,
                     step.increase = 0.024,
                     tip.length = 0.01)+
  theme_bw() +
  ylab("Intensity Mono/Poly ADP-ribose")+
  xlab("Conditions")+
  ggtitle("Zero Mean and Unit Variance Normalization")+
  labs(fill='Mean y-value')

# For easier viewership, and no requriement to explain the negative values, you can add a constant - typically equal to the minimum value.
p_normalised_constant <- ggplot(normalised_combined_df, aes(x=Cell, y = newIntensity488 + abs(min(newIntensity488))))+
  geom_boxplot(aes(fill = mean_newIntensity488 + abs(min(newIntensity488))), outlier.shape = NA) +
  scale_fill_viridis_c(direction = 1)+
  coord_cartesian(ylim = c(0, 8)) +
  stat_pvalue_manual(stat_normalised_NAD,  
                     label = "p.signif", 
                     y.position = 5,
                     step.increase = 0.024,
                     tip.length = 0.01)+
  theme_bw()+
  ylab("Intensity Mono/Poly ADP-ribose")+
  xlab("Conditions")+
  ggtitle("Zero Mean and Unit Variance Normalization (+Constant)")+
  labs(fill='Mean y-value')
p_normalised_constant

