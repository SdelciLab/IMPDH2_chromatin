# libraries
library(tidyverse)
library(ggpubr)

# get data MDAMB231
data_MDA_0h <- read.delim("MDA_NEW0,2,4,24H/0h.imp.parpdoseresp__2023-03-01T07_21_49-Measurement 1/Evaluation8/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

data_MDA_2h <- read.delim("MDA_NEW0,2,4,24H/2h.imp.parp__2023-03-01T09_53_17-Measurement 1/Evaluation5/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

data_MDA_4h <- read.delim("MDA_NEW0,2,4,24H/4h.imp.parp__2023-03-02T07_55_21-Measurement 1/Evaluation3/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

data_MDA_24h <- read.delim("MDA_NEW0,2,4,24H/24.imp.parp__2023-03-02T10_00_39-Measurement 1/Evaluation7/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

# get data Cal51
data_CAL_0h <- read.delim("imp.parpCAL51/T0cal51.imp.parp__2023-03-21T09_12_25-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

data_CAL_4h <- read.delim("imp.parpCAL51/T4cal51.imp.parp1__2023-03-21T07_39_24-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

data_CAL_24h <- read.delim("imp.parpCAL51/T24cal51.imp.parp__2023-03-21T10_53_38-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected Selected.txt", skip = 9)

# explore
colnames(data_MDA_0h)
table(colnames(data_MDA_0h) == colnames(data_MDA_2h))
table(colnames(data_MDA_0h) == colnames(data_MDA_4h))
table(colnames(data_MDA_0h) == colnames(data_MDA_24h))

colnames(data_CAL_0h)
table(colnames(data_CAL_0h) == colnames(data_CAL_4h))
table(colnames(data_CAL_0h) == colnames(data_CAL_24h))

table(data_MDA_0h$Row, data_MDA_0h$Column)
table(data_MDA_2h$Row, data_MDA_2h$Column)
table(data_MDA_4h$Row, data_MDA_4h$Column)
table(data_MDA_24h$Row, data_MDA_24h$Column)

table(data_CAL_0h$Row, data_CAL_0h$Column)
table(data_CAL_4h$Row, data_CAL_4h$Column)
table(data_CAL_24h$Row, data_CAL_24h$Column)

# add time point
data_MDA_0h <- data_MDA_0h %>%
    mutate(Timepoint = rep("0h", nrow(data_MDA_0h)))
data_MDA_2h <- data_MDA_2h %>%
    mutate(Timepoint = rep("2h", nrow(data_MDA_2h)))
data_MDA_4h <- data_MDA_4h %>%
    mutate(Timepoint = rep("4h", nrow(data_MDA_4h)))
data_MDA_24h <- data_MDA_24h %>%
    mutate(Timepoint = rep("24h", nrow(data_MDA_24h)))

data_CAL_0h <- data_CAL_0h %>%
    mutate(Timepoint = rep("0h", nrow(data_CAL_0h)))
data_CAL_4h <- data_CAL_4h %>%
    mutate(Timepoint = rep("4h", nrow(data_CAL_4h)))
data_CAL_24h <- data_CAL_24h %>%
    mutate(Timepoint = rep("24h", nrow(data_CAL_24h)))

# put data together
data_MDA_all <- data_MDA_0h %>%
    rbind(data_MDA_2h) %>%
    rbind(data_MDA_4h) %>%
    rbind(data_MDA_24h)

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

# check
table(data_MDA_all_cond$Condition, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Condition, data_MDA_all_cond$Column)

table(data_MDA_all_cond$Well, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Well, data_MDA_all_cond$Column)

table(data_MDA_all_cond$Timepoint, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Timepoint, data_MDA_all_cond$Column)

# check
table(data_CAL_all_cond$Condition, data_CAL_all_cond$Row)
table(data_CAL_all_cond$Condition, data_CAL_all_cond$Column)

table(data_CAL_all_cond$Well, data_CAL_all_cond$Row)
table(data_CAL_all_cond$Well, data_CAL_all_cond$Column)

table(data_CAL_all_cond$Timepoint, data_CAL_all_cond$Row)
table(data_CAL_all_cond$Timepoint, data_CAL_all_cond$Column)

# factor
data_MDA_all_cond$Timepoint <- factor(data_MDA_all_cond$Timepoint, levels = c("0h", "2h", "4h", "24h"))
data_MDA_all_cond$Condition <- factor(data_MDA_all_cond$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_MDA_all_cond$Well <- factor(data_MDA_all_cond$Well)

data_CAL_all_cond$Timepoint <- factor(data_CAL_all_cond$Timepoint, levels = c("0h", "4h", "24h"))
data_CAL_all_cond$Condition <- factor(data_CAL_all_cond$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_CAL_all_cond$Well <- factor(data_CAL_all_cond$Well)

# normalise IMPDH2 signal by DMSO

# get DMSO every time point
MDA_0h_DMSO <- data_MDA_all_cond %>%
    filter(Timepoint == "0h" & Condition == "DMSO") %>%
    pull(Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean)

MDA_2h_DMSO <- data_MDA_all_cond %>%
    filter(Timepoint == "2h" & Condition == "DMSO") %>%
    pull(Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean)

MDA_4h_DMSO <- data_MDA_all_cond %>%
    filter(Timepoint == "4h" & Condition == "DMSO") %>%
    pull(Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean)

MDA_24h_DMSO <- data_MDA_all_cond %>%
    filter(Timepoint == "24h" & Condition == "DMSO") %>%
    pull(Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean)

# get median
MDA_0h_DMSO_median <- median(MDA_0h_DMSO)
MDA_2h_DMSO_median <- median(MDA_2h_DMSO)
MDA_4h_DMSO_median <- median(MDA_4h_DMSO)
MDA_24h_DMSOmedian <- median(MDA_24h_DMSO)

# get Nuclei.Selected...Intensity.Nucleus.Alexa.488.Mean normalised by DMSO
data_MDA_all_cond_norm <- data_MDA_all_cond %>%
    mutate(IMPDH2_nucl_norm = case_when(
        Timepoint == "0h" ~ Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean/MDA_0h_DMSO_median,
        Timepoint == "2h" ~ Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean/MDA_2h_DMSO_median,
        Timepoint == "4h" ~ Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean/MDA_4h_DMSO_median,
        Timepoint == "24h" ~ Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean/MDA_24h_DMSOmedian)
    )

# check
median(data_MDA_all_cond_norm$IMPDH2_nucl_norm[data_MDA_all_cond_norm$Timepoint == "0h" &
                                                   data_MDA_all_cond_norm$Condition == "DMSO"])
median(data_MDA_all_cond_norm$IMPDH2_nucl_norm[data_MDA_all_cond_norm$Timepoint == "2h" &
                                                   data_MDA_all_cond_norm$Condition == "DMSO"])
median(data_MDA_all_cond_norm$IMPDH2_nucl_norm[data_MDA_all_cond_norm$Timepoint == "4h" &
                                                   data_MDA_all_cond_norm$Condition == "DMSO"])
median(data_MDA_all_cond_norm$IMPDH2_nucl_norm[data_MDA_all_cond_norm$Timepoint == "24h" &
                                                   data_MDA_all_cond_norm$Condition == "DMSO"])

# plot IMPDH2 nuclear in MDAMB231
ggplot(data_MDA_all_cond_norm, aes(x = Timepoint, y = log2(IMPDH2_nucl_norm), fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "Time point", y = "log2(IMPDH2 nuclear mean intensity \n normalised by DMSO)", 
         fill = "") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +                       
    theme_classic()
ggsave("plots/log2IMPDH2nuc_MDA_ETO_normDMSO.pdf", device = "pdf", width = 7, height = 4)

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

p2 <- ggboxplot(
    data_MDA_all_cond_norm, x = "Condition", y = "log2(IMPDH2_nucl_norm)", fill = "Condition",
    facet.by = "Timepoint", nrow = 1) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2
    )

p2 +  labs(y ="log2(IMPDH2 nuclear mean intensity \n normalised by DMSO)") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +  
    scale_y_continuous(limits = c(-2,4.5)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=1))
ggsave("plots/log2IMPDH2nuc_MDA_ETO_norm_DMSOfacetstats.pdf", device = "pdf", width = 7.4, height = 4)

# get number of nuclei in MDAMB231
data_MDA_all_cond_summary <- data_MDA_all_cond %>%
    group_by(Timepoint, Condition, Well) %>%
    summarise(number_of_nuclei = n())

data_MDA_all_cond_summary_stats <- data_MDA_all_cond_summary %>%
    group_by(Timepoint, Condition) %>%
    summarise(n = n(),
              mean_nnuclei = mean(number_of_nuclei),
              sd_nnuclei = sd(number_of_nuclei),
              se_nnuclei = sd_nnuclei/sqrt(n))

# obtain stats 
ttest_pvals_conc <- compare_means(number_of_nuclei ~ Condition, data = data_MDA_all_cond_summary, method = "t.test", alternative = "two.sided", group.by = "Timepoint")
write.csv(ttest_pvals_conc, "t.test_nnuclei_MDAMB231.csv")

# plot
ggplot() +
    geom_col(data = data_MDA_all_cond_summary_stats, aes(x = Condition, y = mean_nnuclei, fill = Condition), 
             width = 0.5) +
    geom_errorbar(data = data_MDA_all_cond_summary_stats, aes(x = Condition, ymin=mean_nnuclei-sd_nnuclei, 
                                                              ymax=mean_nnuclei+sd_nnuclei), colour="black", width=.1) +
    geom_jitter(data = data_MDA_all_cond_summary, aes(x = Condition, y = number_of_nuclei), 
                height = 0, width = 0.1, col = "black")  +
    facet_wrap(~Timepoint) +
    xlab("") +
    ylab("Number of nuclei") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")
ggsave("plots/nnuclei_MDA_etodose.pdf", device = "pdf", width = 5, height =4)

# remove last conc of CAL
data_CAL_all_cond_filt <- data_CAL_all_cond %>%
    filter(Condition != "ETO 10uM")

# get number of nuclei
data_CAL_all_cond_filt_summary <- data_CAL_all_cond_filt %>%
    group_by(Timepoint, Condition, Well) %>%
    summarise(number_of_nuclei = n())

data_CAL_all_cond_filt_summary_stats <- data_CAL_all_cond_filt_summary %>%
    group_by(Timepoint, Condition) %>%
    summarise(n = n(),
              mean_nnuclei = mean(number_of_nuclei),
              sd_nnuclei = sd(number_of_nuclei),
              se_nnuclei = sd_nnuclei/sqrt(n))

# obtain stats 
ttest_pvals_conc <- compare_means(number_of_nuclei ~ Condition, data = data_CAL_all_cond_filt_summary, method = "t.test", alternative = "two.sided", group.by = "Timepoint")
write.csv(ttest_pvals_conc, "t.test_nnuclei_Cal51.csv")

# plot
ggplot() +
    geom_col(data = data_CAL_all_cond_filt_summary_stats, aes(x = Condition, y = mean_nnuclei, fill = Condition), width = 0.5) +
    geom_errorbar(data = data_CAL_all_cond_filt_summary_stats, aes(x = Condition, ymin=mean_nnuclei-sd_nnuclei, ymax=mean_nnuclei+sd_nnuclei), colour="black", width=.1) +
    geom_jitter(data = data_CAL_all_cond_filt_summary, aes(x = Condition, y = number_of_nuclei), height = 0, width = 0.1, col = "black")  +
    facet_wrap(~Timepoint) +
    xlab("") +
    ylab("Number of nuclei") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                                 "ETO 10uM" = "#CC6600")) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")
ggsave("plots/nnuclei_CAL_etodose.pdf", device = "pdf", width = 6, height =3)
