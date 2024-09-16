# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")
library(ggpubr)
library(rstatix)

# get data
data24h_r1_1 <- read.delim("1strep/new24h-spotsrep1/24h.24heto.1__2023-05-29T12_58_47-Measurement 1/Evaluation5/Objects_Population - Cyto Selected.txt", skip = 9)
data24h_spots_r1_1 <- read.delim("1strep/new24h-spotsrep1/24h.24heto.1__2023-05-29T12_58_47-Measurement 1/Evaluation5/Objects_Population - Spots.txt", skip = 9)

data24h_r1_2 <- read.delim("1strep/new24h-spotsrep1/24h.24heto.2__2023-05-29T13_36_40-Measurement 1/Evaluation5/Objects_Population - Cyto Selected.txt", skip = 9)
data24h_spots_r1_2 <- read.delim("1strep/new24h-spotsrep1/24h.24heto.2__2023-05-29T13_36_40-Measurement 1/Evaluation5/Objects_Population - Spots.txt", skip = 9)

data24h_r2 <- read.delim("2ndrep/10uM24h,2ndrepl__2023-06-26T07_11_22-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)
data24h_spots_r2 <- read.delim("2ndrep/10uM24h,2ndrepl__2023-06-26T07_11_22-Measurement 1/Evaluation2/Objects_Population - Spots.txt", skip = 9) 

data24h_r3 <- read.delim("3rdrep/10uM24h,3rdrepl__2023-06-26T09_40_20-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)
data24h_spots_r3 <- read.delim("3rdrep/10uM24h,3rdrepl__2023-06-26T09_40_20-Measurement 1/Evaluation2/Objects_Population - Spots.txt", skip = 9)

# info
colnames(data24h_r1_1)
colnames(data24h_spots_r1_1)

table(data24h_r1_1$Row, data24h_r1_1$Column)
table(data24h_r1_2$Row, data24h_r1_2$Column)

table(data24h_r2$Row, data24h_r2$Column)
table(data24h_r3$Row, data24h_r3$Column)

# check
table(colnames(data24h_r1_1) == colnames(data24h_r1_2))
table(colnames(data24h_r1_1) == colnames(data24h_r2))
table(colnames(data24h_r1_1) == colnames(data24h_r3))

table(colnames(data24h_spots_r1_1) == colnames(data24h_spots_r1_2))
table(colnames(data24h_spots_r1_1) == colnames(data24h_spots_r2))
table(colnames(data24h_spots_r1_1) == colnames(data24h_spots_r3))

# add plate name
data24h_r1_1 <- data24h_r1_1 %>%
    mutate(Replicate = rep(1, nrow(data24h_r1_1)))
data24h_spots_r1_1 <- data24h_spots_r1_1 %>%
    mutate(Replicate = rep(1, nrow(data24h_spots_r1_1)))

data24h_r1_2 <- data24h_r1_2 %>%
    mutate(Replicate = rep(1, nrow(data24h_r1_2)))
data24h_spots_r1_2 <- data24h_spots_r1_2 %>%
    mutate(Replicate = rep(1, nrow(data24h_spots_r1_2)))

data24h_r2 <- data24h_r2 %>%
    mutate(Replicate = rep(2, nrow(data24h_r2)))
data24h_spots_r2 <- data24h_spots_r2 %>%
    mutate(Replicate = rep(2, nrow(data24h_spots_r2)))

data24h_r3 <- data24h_r3 %>%
    mutate(Replicate = rep(3, nrow(data24h_r3)))
data24h_spots_r3 <- data24h_spots_r3 %>%
    mutate(Replicate = rep(3, nrow(data24h_spots_r3)))

# put together data
data24h_all <- data24h_r1_1 %>% 
    rbind(data24h_r1_2) %>% 
    rbind(data24h_r2) %>% 
    rbind(data24h_r3) 

data24h_spots_all <- data24h_spots_r1_1 %>% 
    rbind(data24h_spots_r1_2) %>% 
    rbind(data24h_spots_r2) %>% 
    rbind(data24h_spots_r3) 

# add condition
data24h_all_cond <- data24h_all %>%
    mutate(Cell = case_when(
        Row == 2 | Row == 5 ~ "KOr",
        Row == 3 | Row == 6 ~ "3xnls",
        Row == 4 | Row == 7 ~ "nes"),
        Treatment = case_when(
            (Column == 2 | Column == 3 | Column == 4) & (Row == 2 | Row == 3 | Row == 4)  ~ "DMSO",
            ((Column == 6 | Column == 7 | Column == 8) & (Row == 2 | Row == 3 | Row == 4)) & Replicate != 1  ~ "ETO",
            ((Column == 7 | Column == 8 | Column == 9) & (Row == 2 | Row == 3 | Row == 4)) & Replicate == 1  ~ "ETO",
            (Column == 2 | Column == 3 | Column == 4) & (Row == 5 | Row == 6 | Row == 7)  ~ "MPA",
            ((Column == 6 | Column == 7 | Column == 8) & (Row == 5 | Row == 6 | Row == 7)) & Replicate != 1 ~ "ETO+MPA",
            ((Column == 7 | Column == 8 | Column == 9) & (Row == 5 | Row == 6 | Row == 7)) & Replicate == 1~ "ETO+MPA"),
        Well = case_when(
            (Column == 2 | Column == 6) & Replicate != 1 ~ 1,
            (Column == 3 | Column == 7) & Replicate != 1 ~ 2,
            (Column == 4 | Column == 8) & Replicate != 1 ~ 3,
            (Column == 2 | Column == 7) & Replicate == 1 ~ 1,
            (Column == 3 | Column == 8) & Replicate == 1 ~ 2,
            (Column == 4 | Column == 9) & Replicate == 1 ~ 3)
    )

data24h_spots_all_cond <- data24h_spots_all %>%
    mutate(Cell = case_when(
        Row == 2 | Row == 5 ~ "KOr",
        Row == 3 | Row == 6 ~ "3xnls",
        Row == 4 | Row == 7 ~ "nes"),
        Treatment = case_when(
            (Column == 2 | Column == 3 | Column == 4) & (Row == 2 | Row == 3 | Row == 4)  ~ "DMSO",
            ((Column == 6 | Column == 7 | Column == 8) & (Row == 2 | Row == 3 | Row == 4)) & Replicate != 1  ~ "ETO",
            ((Column == 7 | Column == 8 | Column == 9) & (Row == 2 | Row == 3 | Row == 4)) & Replicate == 1  ~ "ETO",
            (Column == 2 | Column == 3 | Column == 4) & (Row == 5 | Row == 6 | Row == 7)  ~ "MPA",
            ((Column == 6 | Column == 7 | Column == 8) & (Row == 5 | Row == 6 | Row == 7)) & Replicate != 1 ~ "ETO+MPA",
            ((Column == 7 | Column == 8 | Column == 9) & (Row == 5 | Row == 6 | Row == 7)) & Replicate == 1~ "ETO+MPA"),
        Well = case_when(
            (Column == 2 | Column == 6) & Replicate != 1 ~ 1,
            (Column == 3 | Column == 7) & Replicate != 1 ~ 2,
            (Column == 4 | Column == 8) & Replicate != 1 ~ 3,
            (Column == 2 | Column == 7) & Replicate == 1 ~ 1,
            (Column == 3 | Column == 8) & Replicate == 1 ~ 2,
            (Column == 4 | Column == 9) & Replicate == 1 ~ 3)
    )

# check
table(data24h_all_cond$Cell, data24h_all_cond$Row)
table(data24h_all_cond$Cell, data24h_all_cond$Column)
table(data24h_spots_all_cond$Cell, data24h_spots_all_cond$Row)
table(data24h_spots_all_cond$Cell, data24h_spots_all_cond$Column)

table(data24h_all_cond$Treatment, data24h_all_cond$Row)
table(data24h_all_cond$Treatment, data24h_all_cond$Column)
table(data24h_all_cond$Treatment, data24h_all_cond$Row, data24h_all_cond$Replicate)
table(data24h_all_cond$Treatment, data24h_all_cond$Column, data24h_all_cond$Replicate)
table(data24h_spots_all_cond$Treatment, data24h_spots_all_cond$Row)
table(data24h_spots_all_cond$Treatment, data24h_spots_all_cond$Column)

table(data24h_all_cond$Well, data24h_all_cond$Row)
table(data24h_all_cond$Well, data24h_all_cond$Column)
table(data24h_spots_all_cond$Well, data24h_spots_all_cond$Row)
table(data24h_spots_all_cond$Well, data24h_spots_all_cond$Column)

# put together plates

# factors
data24h_all_cond$Cell <- factor(data24h_all_cond$Cell, levels = c("KOr", "3xnls", "nes"))
data24h_all_cond$Treatment <- factor(data24h_all_cond$Treatment, levels = c("DMSO", "ETO", "MPA", "ETO+MPA"))
data24h_all_cond$Well <- factor(data24h_all_cond$Well)
data24h_all_cond$Replicate <- factor(data24h_all_cond$Replicate)

# add logs intensities
data24h_all_cond_log <- data24h_all_cond %>%
    mutate(IMPDH2_cyt_log2 = log2(Cyto.Selected...Intensity.Cytoplasm.Alexa.647.Mean),
           IMPDH2_nuc_log2 = log2(Cyto.Selected...Intensity.Nucleus.Alexa.647..Mean),
           PARP1_cyt_log2 = log2(Cyto.Selected...Intensity.Cytoplasm.Alexa.488.Mean),
           PARP1_nuc_log2 = log2(Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean))

# remove outlier
data24h_all_cond_log_noout <- data24h_all_cond_log %>%
    filter(Cyto.Selected...Formula < 6 & Cyto.Selected...Total.Spot.Area < 3000 & 
               Cyto.Selected...Number.of.Spots < 150 & Cyto.Selected...Spots.Area..µm.. < 900)

# only DMSO ETO
data24h_all_cond_log_noout_DMSO_ETO <- data24h_all_cond_log_noout %>%
    filter(Treatment %in% c("DMSO", "ETO"))

# only DMSO ETO
data24h_all_cond_log_noout_DMSO_ETO <- data24h_all_cond_log_noout %>%
    filter(Treatment %in% c("DMSO", "ETO"))

# plot number of PARP1 spots
ggplot(data24h_all_cond_log_noout_DMSO_ETO, aes(x = Cell, y = log2(Cyto.Selected...Number.of.Spots), fill = Treatment)) +
    geom_boxplot(position = "dodge2", fatten = 2, width = 0.75) +
    labs(x = "", y = "log2(Number of cytosolic PARP1 spots)", fill = "") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC","ETO" = "#CC6600")) +  
    theme_classic() +
    scale_x_discrete(labels=c("KOr" = "KO+WT", "3xnls" = "KO+NLS",
                              "nes" = "KO+NES")) +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.text = element_text(size = 11), axis.title = element_text(size=12))
ggsave("plots/final_24h_DMSO_ETO_nspots_boxplot.pdf", device = "pdf", width = 5, height = 3.5)

# facet
ggplot(data24h_all_cond_log_noout_DMSO_ETO, aes(x = Treatment, y = log2(Cyto.Selected...Number.of.Spots), fill = Treatment)) +
    facet_wrap(~Cell, labeller = labeller(Cell = c("KOr" = "KO+WT", "3xnls" = "KO+NLS",
                                                   "nes" = "KO+NES"))) +
    geom_boxplot(position = "dodge2", fatten = 2, width = 0.75) +
    labs(x = "", y = "Number of cytosolic PARP1 spots", fill = "") +
    scale_fill_manual(name = "", labels = c("DMSO", "ETO"),
                      values = c("DMSO" = "#FFE5CC","ETO" = "#CC6600"))+
    stat_compare_means(method = "wilcox", comparisons = list(c("DMSO","ETO")), size = 2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave("plots/final_24h_DMSO_ETO_nspots_facet.pdf", device = "pdf", width = 5, height = 3)

# by treatments
ggplot(data24h_all_cond_log_noout_DMSO_ETO, aes(x = Cell, y = log2(Cyto.Selected...Number.of.Spots), fill = Cell)) +
    facet_wrap(~Treatment, labeller = labeller(Cell = c("KOr" = "KO+WT", "3xnls" = "KO+NLS",
                                                        "nes" = "KO+NES"))) +
    geom_boxplot(position = "dodge2", fatten = 2, width = 0.75) +
    labs(x = "", y = "Number of cytosolic PARP1 spots", fill = "") +
    scale_fill_manual(name = "", labels = c("KO-WT", "KO-NLS", "KO-NES"),
                      values = c("KOr" = "#4C0099", "3xnls" = "#855FFF", "nes" = "#A392C0"))+
    stat_compare_means(method = "wilcox", comparisons = list(c("KOr","3xnls"), c("3xnls", "nes"), c("KOr", "nes")), size = 2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave("plots/final_24h_cells_nspots_facet.pdf", device = "pdf", width = 5, height = 3)

# numbers
table(data24h_all_cond_log_noout_DMSO_ETO$Treatment, data24h_all_cond_log_noout_DMSO_ETO$Cell)
