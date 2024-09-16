# libraries
library(tidyverse)
library(ggpubr)

# get data MDAMB231
data_MDA_0h_1 <- read.delim("MDA-231/MDAWT0H-1__2023-10-27T12_38_08-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_0h_2 <- read.delim("MDA-231/MDAWT0H-2__2023-10-27T13_34_04-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_2h_1 <- read.delim("MDA-231/MDAWT2H-1__2023-10-27T11_30_29-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_2h_2 <- read.delim("MDA-231/MDAWT2H-2__2023-10-27T12_09_12-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_4h_1 <- read.delim("MDA-231/MDAWT4H-1__2023-10-27T10_50_05-Measurement 1/Evaluation3/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_4h_2 <- read.delim("MDA-231/MDAWT4H-2__2023-10-27T11_06_20-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_24h_1 <- read.delim("MDA-231/MDAWT24H-1__2023-10-27T10_13_03-Measurement 2/Evaluation1/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_24h_2 <- read.delim("MDA-231/MDAWT24H-2__2023-10-27T10_33_03-Measurement 1/Evaluation1/Objects_Population - Cyto Selected.txt", skip = 9)

# get data Cal51
data_CAL_0h <- read.delim("cal51.1/T0cal51.imp.parp__2023-03-21T09_12_25-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_CAL_4h <- read.delim("cal51.1/T4cal51.imp.parp1__2023-03-21T07_39_24-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)

data_CAL_24h <- read.delim("cal51.1/T24cal51.imp.parp__2023-03-21T10_53_38-Measurement 1/Evaluation3/Objects_Population - Cyto Selected.txt", skip = 9)

# explore
colnames(data_MDA_0h_1)
table(colnames(data_MDA_0h_1) == colnames(data_MDA_0h_2))
table(colnames(data_MDA_0h_1) == colnames(data_MDA_2h_1))
table(colnames(data_MDA_0h_1) == colnames(data_MDA_2h_2))
table(colnames(data_MDA_0h_1) == colnames(data_MDA_4h_1))
table(colnames(data_MDA_0h_1) == colnames(data_MDA_4h_2))
table(colnames(data_MDA_0h_1) == colnames(data_MDA_24h_1))
table(colnames(data_MDA_0h_1) == colnames(data_MDA_24h_2))

colnames(data_CAL_0h)
table(colnames(data_CAL_0h) == colnames(data_CAL_4h))
table(colnames(data_CAL_0h) == colnames(data_CAL_24h))

table(data_MDA_0h_1$Row, data_MDA_0h_1$Column)
table(data_MDA_0h_2$Row, data_MDA_0h_2$Column)
table(data_MDA_2h_1$Row, data_MDA_2h_1$Column)
table(data_MDA_2h_2$Row, data_MDA_2h_2$Column)
table(data_MDA_4h_1$Row, data_MDA_4h_1$Column)
table(data_MDA_4h_2$Row, data_MDA_4h_2$Column)
table(data_MDA_24h_1$Row, data_MDA_24h_1$Column)
table(data_MDA_24h_2$Row, data_MDA_24h_2$Column)

table(data_CAL_0h$Row, data_CAL_0h$Column)
table(data_CAL_4h$Row, data_CAL_4h$Column)
table(data_CAL_24h$Row, data_CAL_24h$Column)

# add time point and plate
data_MDA_0h_1 <- data_MDA_0h_1 %>%
    mutate(Timepoint = rep("0h", nrow(data_MDA_0h_1)),
           Plate = rep(1, nrow(data_MDA_0h_1)))
data_MDA_0h_2 <- data_MDA_0h_2 %>%
    mutate(Timepoint = rep("0h", nrow(data_MDA_0h_2)),
           Plate = rep(2, nrow(data_MDA_0h_2)))
data_MDA_2h_1 <- data_MDA_2h_1 %>%
    mutate(Timepoint = rep("2h", nrow(data_MDA_2h_1)),
           Plate = rep(1, nrow(data_MDA_2h_1)))
data_MDA_2h_2 <- data_MDA_2h_2 %>%
    mutate(Timepoint = rep("2h", nrow(data_MDA_2h_2)),
           Plate = rep(2, nrow(data_MDA_2h_2)))
data_MDA_4h_1 <- data_MDA_4h_1 %>%
    mutate(Timepoint = rep("4h", nrow(data_MDA_4h_1)),
           Plate = rep(1, nrow(data_MDA_4h_1)))
data_MDA_4h_2 <- data_MDA_4h_2 %>%
    mutate(Timepoint = rep("4h", nrow(data_MDA_4h_2)),
           Plate = rep(2, nrow(data_MDA_4h_2)))
data_MDA_24h_1 <- data_MDA_24h_1 %>%
    mutate(Timepoint = rep("24h", nrow(data_MDA_24h_1)),
           Plate = rep(1, nrow(data_MDA_24h_1)))
data_MDA_24h_2 <- data_MDA_24h_2 %>%
    mutate(Timepoint = rep("24h", nrow(data_MDA_24h_2)),
           Plate = rep(2, nrow(data_MDA_24h_2)))

data_CAL_0h <- data_CAL_0h %>%
    mutate(Timepoint = rep("0h", nrow(data_CAL_0h)))
data_CAL_4h <- data_CAL_4h %>%
    mutate(Timepoint = rep("4h", nrow(data_CAL_4h)))
data_CAL_24h <- data_CAL_24h %>%
    mutate(Timepoint = rep("24h", nrow(data_CAL_24h)))

# put data together
data_MDA_all <- data_MDA_0h_1 %>%
    rbind(data_MDA_0h_2) %>%
    rbind(data_MDA_2h_1) %>%
    rbind(data_MDA_2h_2) %>%
    rbind(data_MDA_4h_1) %>%
    rbind(data_MDA_4h_2) %>%
    rbind(data_MDA_24h_1) %>%
    rbind(data_MDA_24h_2)

data_CAL_all <- data_CAL_0h %>%
    rbind(data_CAL_4h) %>%
    rbind(data_CAL_24h) 

# add conditions
data_MDA_all_cond <- data_MDA_all %>%
    mutate(Condition = case_when(
        Column == 3 ~ "DMSO",
        Column == 4 ~ "ETO 1uM",
        Column == 5 ~ "ETO 2.5uM",
        Column == 6 ~ "ETO 5uM",
        Column == 7 ~ "ETO 10uM"),
        Well = case_when(
            Row == 5 ~ 1, 
            Row == 6 ~ 2, 
            Row == 7 ~ 3)
    )

data_CAL_all_cond <- data_CAL_all %>%
    mutate(Condition = case_when(
        Column == 3 ~ "DMSO",
        Column == 4 ~ "ETO 1uM",
        Column == 5 ~ "ETO 2.5uM",
        Column == 6 ~ "ETO 5uM",
        Column == 7 ~ "ETO 10uM"),
        Well = case_when(
            Row == 2 ~ 1, 
            Row == 3 ~ 2, 
            Row == 4 ~ 3)
    )

# remove one well: G5 from MDA time 0
data_MDA_all_cond2 <- data_MDA_all_cond %>%
    filter(!(Timepoint == "0h" & Condition == "ETO 2.5uM" &  Well == 3))

# check
table(data_MDA_all_cond2$Condition, data_MDA_all_cond2$Row)
table(data_MDA_all_cond2$Condition, data_MDA_all_cond2$Column)

table(data_MDA_all_cond2$Well, data_MDA_all_cond2$Row)
table(data_MDA_all_cond2$Well, data_MDA_all_cond2$Column)

table(data_MDA_all_cond2$Timepoint, data_MDA_all_cond2$Row)
table(data_MDA_all_cond2$Timepoint, data_MDA_all_cond2$Column)

table(data_MDA_all_cond2$Plate, data_MDA_all_cond2$Row)
table(data_MDA_all_cond2$Plate, data_MDA_all_cond2$Column)

table(data_MDA_all_cond2$Timepoint, data_MDA_all_cond2$Condition, 
      data_MDA_all_cond2$Well, data_MDA_all_cond2$Plate)

# check
table(data_CAL_all_cond$Condition, data_CAL_all_cond$Row)
table(data_CAL_all_cond$Condition, data_CAL_all_cond$Column)

table(data_CAL_all_cond$Well, data_CAL_all_cond$Row)
table(data_CAL_all_cond$Well, data_CAL_all_cond$Column)

table(data_CAL_all_cond$Timepoint, data_CAL_all_cond$Row)
table(data_CAL_all_cond$Timepoint, data_CAL_all_cond$Column)

# factor
data_MDA_all_cond2$Timepoint <- factor(data_MDA_all_cond2$Timepoint, levels = c("0h", "2h", "4h", "24h"))
data_MDA_all_cond2$Condition <- factor(data_MDA_all_cond2$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_MDA_all_cond2$Well <- factor(data_MDA_all_cond2$Well)
data_MDA_all_cond2$Plate <- factor(data_MDA_all_cond2$Plate)

data_CAL_all_cond$Timepoint <- factor(data_CAL_all_cond$Timepoint, levels = c("0h", "4h", "24h"))
data_CAL_all_cond$Condition <- factor(data_CAL_all_cond$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_CAL_all_cond$Well <- factor(data_CAL_all_cond$Well)

# MDAMB231

# normalise number of PARP1 spots by DMSO

# get DMSO every time point
MDA_0h_DMSO <- data_MDA_all_cond2 %>%
    filter(Timepoint == "0h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Number.of.Spots)

MDA_2h_DMSO <- data_MDA_all_cond2 %>%
    filter(Timepoint == "2h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Number.of.Spots)

MDA_4h_DMSO <- data_MDA_all_cond2 %>%
    filter(Timepoint == "4h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Number.of.Spots)

MDA_24h_DMSO <- data_MDA_all_cond2 %>%
    filter(Timepoint == "24h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Number.of.Spots)

# get median
MDA_0h_DMSO_median <- median(MDA_0h_DMSO, na.rm = T)
MDA_2h_DMSO_median <- median(MDA_2h_DMSO, na.rm = T)
MDA_4h_DMSO_median <- median(MDA_4h_DMSO, na.rm = T)
MDA_24h_DMSOmedian <- median(MDA_24h_DMSO, na.rm = T)

# get Cyto.Selected...Number.of.Spots normalised by DMSO
data_MDA_all_cond_norm <- data_MDA_all_cond2 %>%
    mutate(PARP1_nspots_norm = case_when(
        Timepoint == "0h" ~ Cyto.Selected...Number.of.Spots/MDA_0h_DMSO_median,
        Timepoint == "2h" ~ Cyto.Selected...Number.of.Spots/MDA_2h_DMSO_median,
        Timepoint == "4h" ~ Cyto.Selected...Number.of.Spots/MDA_4h_DMSO_median,
        Timepoint == "24h" ~ Cyto.Selected...Number.of.Spots/MDA_24h_DMSOmedian)
    )

# check
median(data_MDA_all_cond_norm$PARP1_nspots_norm[data_MDA_all_cond_norm$Timepoint == "0h" &
                                                    data_MDA_all_cond_norm$Condition == "DMSO"], na.rm = T)
median(data_MDA_all_cond_norm$PARP1_nspots_norm[data_MDA_all_cond_norm$Timepoint == "2h" &
                                                    data_MDA_all_cond_norm$Condition == "DMSO"], na.rm = T)
median(data_MDA_all_cond_norm$PARP1_nspots_norm[data_MDA_all_cond_norm$Timepoint == "4h" &
                                                    data_MDA_all_cond_norm$Condition == "DMSO"], na.rm = T)
median(data_MDA_all_cond_norm$PARP1_nspots_norm[data_MDA_all_cond_norm$Timepoint == "24h" &
                                                    data_MDA_all_cond_norm$Condition == "DMSO"], na.rm = T)

# plot number of PARP1 spots normalised by DMSO
ggplot(data_MDA_all_cond_norm, aes(x = Timepoint, y = PARP1_nspots_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "Time point", y = "Number of PARP1 spots \n normalised by DMSO)", 
         fill = "") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +                       
    theme_classic()
ggsave("plots/PARP1_nspots_MDA_ETO_normDMSO.pdf", device = "pdf", width = 7, height = 4)

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

PARP1nspotsnorm <- ggboxplot(
    data_MDA_all_cond_norm, x = "Condition", y = "PARP1_nspots_norm", fill = "Condition",
    facet.by = "Timepoint", nrow = 1) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2
    )

PARP1nspotsnorm +  labs(y ="Number of PARP1 spots \n normalised by DMSO)") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=1))
ggsave("plots/PARP1nspots_MDA_ETO_norm_DMSOfacetstats.pdf", device = "pdf", width = 7.4, height = 4)

# CAL51

# remove last conc of CAL
data_CAL_all_cond_filt <- data_CAL_all_cond %>%
    filter(Condition != "ETO 10uM")

# normalise IMPDH2 intensity by DMSO condition

# get DMSO every time point
CAL_0h_DMSO <- data_CAL_all_cond %>%
    filter(Timepoint == "0h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean)

CAL_4h_DMSO <- data_CAL_all_cond %>%
    filter(Timepoint == "4h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean)

CAL_24h_DMSO <- data_CAL_all_cond %>%
    filter(Timepoint == "24h" & Condition == "DMSO") %>%
    pull(Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean)

# get median
CAL_0h_DMSO_median <- median(CAL_0h_DMSO)
CAL_4h_DMSO_median <- median(CAL_4h_DMSO)
CAL_24h_DMSO_median <- median(CAL_24h_DMSO)

# get Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean normalised by DMSO
data_CAL_all_cond_norm <- data_CAL_all_cond %>%
    mutate(IMPDH2_nucl_norm = case_when(
        Timepoint == "0h" ~ Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean/CAL_0h_DMSO_median,
        Timepoint == "4h" ~ Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean/CAL_4h_DMSO_median,
        Timepoint == "24h" ~ Cyto.Selected...Intensity.Nucleus.Alexa.488.Mean/CAL_24h_DMSO_median)
    )

# check
median(data_CAL_all_cond_norm$IMPDH2_nucl_norm[data_CAL_all_cond_norm$Timepoint == "0h" &
                                                   data_CAL_all_cond_norm$Condition == "DMSO"])
median(data_CAL_all_cond_norm$IMPDH2_nucl_norm[data_CAL_all_cond_norm$Timepoint == "4h" &
                                                   data_CAL_all_cond_norm$Condition == "DMSO"])
median(data_CAL_all_cond_norm$IMPDH2_nucl_norm[data_CAL_all_cond_norm$Timepoint == "24h" &
                                                   data_CAL_all_cond_norm$Condition == "DMSO"])

# plot IMPDH2 nuclear
ggplot(data_CAL_all_cond_norm, aes(x = Timepoint, y = log2(IMPDH2_nucl_norm), fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "Time point", y = "log2(IMPDH2 nuclear mean intensity \n normalised by DMSO)", 
         fill = "") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +                       
    theme_classic()
ggsave("plots/log2IMPDH2nuc_CAL_ETO_normDMSO.pdf", device = "pdf", width = 7, height = 4)

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

pcal51_IMPDH_norm <- ggboxplot(
    data_CAL_all_cond_norm, x = "Condition", y = "log2(IMPDH2_nucl_norm)", fill = "Condition",
    facet.by = "Timepoint", nrow = 1) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2
    )

pcal51_IMPDH_norm +  labs(y ="log2(IMPDH2 nuclear mean intensity \n normalised by DMSO)") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=1))
ggsave("plots/log2IMPDH2nuc_CAL_ETO_norm_DMSOfacetstats.pdf", device = "pdf", width = 7.4, height = 4)
