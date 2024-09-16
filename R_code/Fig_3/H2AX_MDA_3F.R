# libraries
library(tidyverse)
library(ggpubr)

# get data MDAMB231
data_MDA_0h <- read.delim("FOCIgH2AX.impdh2/0.imp.parpTG__2023-03-03T13_17_05-Measurement 1/Evaluation1/Objects_Population - full nuclei.txt", skip = 9)

data_MDA_2h <- read.delim("FOCIgH2AX.impdh2/2.imp.parpTG__2023-03-03T12_35_23-Measurement 1/Evaluation1/Objects_Population - full nuclei.txt", skip = 9)

data_MDA_4h <- read.delim("FOCIgH2AX.impdh2/4.imp.parpTG__2023-03-03T09_44_15-Measurement 1/Evaluation1/Objects_Population - full nuclei.txt", skip = 9)

data_MDA_24h <- read.delim("FOCIgH2AX.impdh2/24.imp.parpTG__2023-03-03T09_09_12-Measurement 1/Evaluation1/Objects_Population - full nuclei.txt", skip = 9)

# explore
colnames(data_MDA_0h)
table(colnames(data_MDA_0h) == colnames(data_MDA_2h))
table(colnames(data_MDA_0h) == colnames(data_MDA_4h))
table(colnames(data_MDA_0h) == colnames(data_MDA_24h))

table(data_MDA_0h$Row, data_MDA_0h$Column)
table(data_MDA_2h$Row, data_MDA_2h$Column)
table(data_MDA_4h$Row, data_MDA_4h$Column)
table(data_MDA_24h$Row, data_MDA_24h$Column)

# add time point
data_MDA_0h <- data_MDA_0h %>%
    mutate(Timepoint = rep("0h", nrow(data_MDA_0h)))
data_MDA_2h <- data_MDA_2h %>%
    mutate(Timepoint = rep("2h", nrow(data_MDA_2h)))
data_MDA_4h <- data_MDA_4h %>%
    mutate(Timepoint = rep("4h", nrow(data_MDA_4h)))
data_MDA_24h <- data_MDA_24h %>%
    mutate(Timepoint = rep("24h", nrow(data_MDA_24h)))

# put data together
data_MDA_all <- data_MDA_0h %>%
    rbind(data_MDA_2h) %>%
    rbind(data_MDA_4h) %>%
    rbind(data_MDA_24h)

# add conditions
data_MDA_all_cond <- data_MDA_all %>%
    mutate(Condition = case_when(
        Column == 3 ~ "DMSO",
        Column == 4 ~ "ETO 1uM",
        Column == 5 ~ "ETO 2.5uM",
        Column == 6 ~ "ETO 5uM",
        Column == 7 ~ "ETO 10uM"),
        Well = case_when(
            Row == 2 ~ 1, 
            Row == 3 ~ 2, 
            Row == 4 ~ 3),
        n_spots_mod = full.nuclei...Number.of.Spots + 0.01
    )

# check
table(data_MDA_all_cond$Condition, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Condition, data_MDA_all_cond$Column)

table(data_MDA_all_cond$Well, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Well, data_MDA_all_cond$Column)

table(data_MDA_all_cond$Timepoint, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Timepoint, data_MDA_all_cond$Column)

# factor
data_MDA_all_cond$Timepoint <- factor(data_MDA_all_cond$Timepoint, levels = c("0h", "2h", "4h", "24h"))
data_MDA_all_cond$Condition <- factor(data_MDA_all_cond$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_MDA_all_cond$Well <- factor(data_MDA_all_cond$Well)

# filter spots higher than 100
data_MDA_all_cond_filt <- data_MDA_all_cond %>%
    filter(full.nuclei...Number.of.Spots < 100)

# normalise H2AX spots by DMSO condition

# get DMSO every time point
MDA_0h_DMSO <- data_MDA_all_cond_filt %>%
    filter(Timepoint == "0h" & Condition == "DMSO") %>%
    pull(full.nuclei...Number.of.Spots)

MDA_2h_DMSO <- data_MDA_all_cond_filt %>%
    filter(Timepoint == "2h" & Condition == "DMSO") %>%
    pull(full.nuclei...Number.of.Spots)

MDA_4h_DMSO <- data_MDA_all_cond_filt %>%
    filter(Timepoint == "4h" & Condition == "DMSO") %>%
    pull(full.nuclei...Number.of.Spots)

MDA_24h_DMSO <- data_MDA_all_cond_filt %>%
    filter(Timepoint == "24h" & Condition == "DMSO") %>%
    pull(full.nuclei...Number.of.Spots)

# get median 
MDA_0h_DMSO_median <- median(MDA_0h_DMSO)
MDA_2h_DMSO_median <- median(MDA_2h_DMSO)
MDA_4h_DMSO_median <- median(MDA_4h_DMSO)
MDA_24h_DMSO_median <- median(MDA_24h_DMSO)

# get full.nuclei...Number.of.Spots normalised by DMSO
data_MDA_all_cond_norm <- data_MDA_all_cond_filt %>%
    mutate(H2AXspots_norm = case_when(
        Timepoint == "0h" ~ full.nuclei...Number.of.Spots/MDA_0h_DMSO_median,
        Timepoint == "2h" ~ full.nuclei...Number.of.Spots/MDA_2h_DMSO_median,
        Timepoint == "4h" ~ full.nuclei...Number.of.Spots/MDA_4h_DMSO_median,
        Timepoint == "24h" ~ full.nuclei...Number.of.Spots/MDA_24h_DMSO_median)
    )

# check
median(data_MDA_all_cond_norm$H2AXspots_norm[data_MDA_all_cond_norm$Timepoint == "0h" &
                                                 data_MDA_all_cond_norm$Condition == "DMSO"])
median(data_MDA_all_cond_norm$H2AXspots_norm[data_MDA_all_cond_norm$Timepoint == "2h" &
                                                 data_MDA_all_cond_norm$Condition == "DMSO"])
median(data_MDA_all_cond_norm$H2AXspots_norm[data_MDA_all_cond_norm$Timepoint == "4h" &
                                                 data_MDA_all_cond_norm$Condition == "DMSO"])
median(data_MDA_all_cond_norm$H2AXspots_norm[data_MDA_all_cond_norm$Timepoint == "24h" &
                                                 data_MDA_all_cond_norm$Condition == "DMSO"])

# plot H2AX spots
ggplot(data_MDA_all_cond_norm, aes(x = Timepoint, y = H2AXspots_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "Time point", y = "Number of H2AX foci \n normalised by DMSO", 
         fill = "") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +                       
    theme_classic()
ggsave("plots/H2AXnspots_MDA_ETO_normDMSO.pdf", device = "pdf", width = 7, height = 4)

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

p2 <- ggboxplot(
    data_MDA_all_cond_norm, x = "Condition", y = "H2AXspots_norm", fill = "Condition",
    facet.by = "Timepoint", nrow = 1) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2
    )

p2 +  labs(y ="Number of H2AX foci \n normalised by DMSO") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=1))
ggsave("plots/H2AXnspots_MDA_ETO_norm_DMSOfacetstats.pdf", device = "pdf", width = 7.4, height = 4)