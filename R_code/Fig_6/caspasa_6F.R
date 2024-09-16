# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")
library(ggpubr)
library(gtools)

# get data
data_casp_r1 <- read.delim("CLEAVEDCASP3/Recconst.cleavedcasp3 rep1__2023-06-27T07_36_42-Measurement 1/Evaluation7/Objects_Population - Cyto Selected.txt", skip = 9)

data_casp_r2 <- read.delim("CLEAVEDCASP3/Recconst.cleavedcasp3 rep2__2023-06-27T11_40_49-Measurement 1/Evaluation1/Objects_Population - Cyto Selected.txt", skip = 9)

# info
colnames(data_casp_r1)

table(data_casp_r1$Row, data_casp_r1$Column)
table(data_casp_r2$Row, data_casp_r2$Column)

# check
table(colnames(data_casp_r1) == colnames(data_casp_r2))

# add plate name 
data_casp_r1 <- data_casp_r1 %>%
    mutate(Plate = rep(1, nrow(data_casp_r1)))
data_casp_r2 <- data_casp_r2 %>%
    mutate(Plate = rep(2, nrow(data_casp_r2)))

# put together data
data_casp_all <- data_casp_r1 %>% 
    rbind(data_casp_r2) 

# add condition
data_casp_all_cond <- data_casp_all %>%
    mutate(Cell = case_when(
        (Row == 2 | Row == 3) &
            (Column == 2 | Column == 3 | Column == 4 | Column == 5 | Column == 6) ~ "KOr",
        (Row == 2 | Row == 3) &
            (Column == 7 | Column == 8 | Column == 9 | Column == 10 | Column == 11) ~ "3xnls",
        (Row == 5 | Row == 6) &
            (Column == 4 | Column == 5 | Column == 6 | Column == 7 | Column == 8) ~ "nes"))

data_casp_all_cond_cond2 <- data_casp_all_cond %>%
    mutate(Treatment = case_when(
        (Cell == "KOr") & (Column == 2 | Column == 3) ~ "DMSO",
        (Cell == "3xnls") & (Column == 7 | Column == 8) ~ "DMSO",
        (Cell == "nes") & (Column == 4 | Column == 5) ~ "DMSO",
        (Cell == "KOr") & (Column == 4 | Column == 5 | Column == 6) ~ "ETO",
        (Cell == "3xnls") & (Column == 9 | Column == 10 | Column == 11) ~ "ETO",
        (Cell == "nes") & (Column == 6 | Column == 7 | Column == 8) ~ "ETO"),
        Well = case_when(
            (Cell == "KOr") & (Column == 2 | Column == 4) ~ 1,
            (Cell == "3xnls") & (Column == 7 | Column == 9) ~ 1,
            (Cell == "nes") & (Column == 4 | Column == 6) ~ 1,
            (Cell == "KOr") & (Column == 3 | Column == 5) ~ 2,
            (Cell == "3xnls") & (Column == 8 | Column == 10) ~ 2,
            (Cell == "nes") & (Column == 5 | Column == 7) ~ 2,
            (Cell == "KOr") & (Column == 6) ~ 3,
            (Cell == "3xnls") & (Column == 11) ~ 3,
            (Cell == "nes") & (Column == 8) ~ 3)
    )

# check
table(data_casp_all_cond_cond2$Cell, data_casp_all_cond_cond2$Row)
table(data_casp_all_cond_cond2$Cell, data_casp_all_cond_cond2$Column)

table(data_casp_all_cond_cond2$Treatment, data_casp_all_cond_cond2$Row)
table(data_casp_all_cond_cond2$Treatment, data_casp_all_cond_cond2$Column)

table(data_casp_all_cond_cond2$Well, data_casp_all_cond_cond2$Row)
table(data_casp_all_cond_cond2$Well, data_casp_all_cond_cond2$Column)

# factors
data_casp_all_cond_cond2$Cell <- factor(data_casp_all_cond_cond2$Cell, levels = c("KOr", "3xnls", "nes"))
data_casp_all_cond_cond2$Treatment <- factor(data_casp_all_cond_cond2$Treatment, levels = c("DMSO", "ETO"))
data_casp_all_cond_cond2$Well <- factor(data_casp_all_cond_cond2$Well)
data_casp_all_cond_cond2$Plate <- factor(data_casp_all_cond_cond2$Plate)

# remove outliers
data_casp_all_cond_cond2_nout <- data_casp_all_cond_cond2 %>%
    filter(Cyto.Selected...Total.Spot.Area < 4000 &
               Cyto.Selected...Spots.Area..µm.. < 900 & 
               Cyto.Selected...Formula < 6)

# check
table(data_casp_all_cond_cond2$Cyto.Selected...Total.Spot.Area > 4000 | 
          data_casp_all_cond_cond2$Cyto.Selected...Spots.Area..µm.. > 900 |
          data_casp_all_cond_cond2$Cyto.Selected...Formula > 6)

# consider wells as replicates
data_casp_all_cond_reps <- data_casp_all_cond_cond2_nout %>%
    mutate(Replicate = case_when(
        (Plate == 1 & Treatment == "DMSO" & Well==1) ~ 1,
        (Plate == 1 & Treatment == "DMSO" & Well==2) ~ 2,
        (Plate == 2 & Treatment == "DMSO" & Well==1) ~ 3,
        (Plate == 2 & Treatment == "DMSO" & Well==2) ~ 4,
        (Plate == 1 & Treatment == "ETO" & Well==1) ~ 1,
        (Plate == 1 & Treatment == "ETO" & Well==2) ~ 2,
        (Plate == 1 & Treatment == "ETO" & Well==3) ~ 3,
        (Plate == 2 & Treatment == "ETO" & Well==1) ~ 4,
        (Plate == 2 & Treatment == "ETO" & Well==2) ~ 5,
        (Plate == 2 & Treatment == "ETO" & Well==3) ~ 6)
    )

# check
table(data_casp_all_cond_reps$Replicate, data_casp_all_cond_reps$Well,
      data_casp_all_cond_reps$Plate, data_casp_all_cond_reps$Treatment)

# plot cytosolic total spot area
ggplot(data_casp_all_cond_reps, aes(x = Cell, y = log2(Cyto.Selected...Total.Spot.Area), fill = Treatment)) +
    geom_boxplot(position = "dodge2", fatten = 2) +
    labs(x = "", y = "log2(Total Cytosolic Casp3 Spot Area)", fill = "") +
    scale_fill_manual(name = "", labels = c("DMSO", "ETO"),
                      values = c("DMSO" = "#FFE5CC","ETO" = "#CC6600")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave("plots/final_total_spot_area_treats_boxplot.pdf", device = "pdf", width = 5, height = 3)

# facet stats
ggplot(data_casp_all_cond_reps, aes(x = Treatment, y = log2(Cyto.Selected...Total.Spot.Area), fill = Treatment)) +
    facet_wrap(~Cell, labeller = labeller(Cell = c("KOr" = "KO+WT", "3xnls" = "KO+NLS",
                                                   "nes" = "KO+NES"))) +
    geom_boxplot(position = "dodge2", fatten = 2, width = 0.75) +
    labs(x = "", y = "log2(Total Cytosolic Casp3 Spot Area)", fill = "") +
    scale_fill_manual(name = "", labels = c("DMSO", "ETO"),
                      values = c("DMSO" = "#FFE5CC","ETO" = "#CC6600"))+
    stat_compare_means(method = "wilcox", comparisons = list(c("DMSO","ETO")), size = 2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave("plots/final_total_spot_area_treats_boxplot_facet.pdf", device = "pdf", width = 5, height = 3)

ggplot(data_casp_all_cond_reps, aes(x = Cell, y = log2(Cyto.Selected...Total.Spot.Area), fill = Cell)) +
    facet_wrap(~Treatment, labeller = labeller(Cell = c("KOr" = "KO+WT", "3xnls" = "KO+NLS",
                                                        "nes" = "KO+NES"))) +
    geom_boxplot(position = "dodge2", fatten = 2, width = 0.75) +
    labs(x = "", y = "log2(Total Cytosolic Casp3 Spot Area)", fill = "") +
    scale_fill_manual(name = "", labels = c("KO-WT", "KO-NLS", "KO-NES"),
                      values = c("KOr" = "#4C0099", "3xnls" = "#855FFF", "nes" = "#A392C0"))+
    stat_compare_means(method = "wilcox", comparisons = list(c("KOr","3xnls"), c("3xnls", "nes"), c("KOr", "nes")), size = 2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave("plots/final_total_spot_area_cells_boxplot_facet.pdf", device = "pdf", width = 5, height = 3)

# numbers
table(data_casp_all_cond_reps$Cell, data_casp_all_cond_reps$Treatment)