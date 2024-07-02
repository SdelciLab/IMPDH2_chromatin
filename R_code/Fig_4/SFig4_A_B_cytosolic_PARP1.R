# Conditions
#2 cell lines:
#* MDA-MB-231
#* Cal51

#Cell lines in different plates.
#Time-points in different plates.

#Scheme of MDA-MB-231:
  
#Treatment:
  
#* DMSO
#* ETO 1 uM
#* ETO 2.5 uM
#* ETO 5 uM
#* ETO 10 uM

#Time-points:
  
#* 0h
#* 2h
#* 4h
#* 24h

#Treatment in columns:
  
#* DMSO: 3
#* ETO 1 uM: 4
#* ETO 2.5 uM: 5
#* ETO 5 uM: 6
#* ETO 10 uM: 7

#Wells in rows:
  
#* 1: 5
#* 2: 6
#* 3: 7

#Scheme of Cal51:
  
 # Treatment:
  
  #* DMSO
#* ETO 1 uM
#* ETO 2.5 uM
#* ETO 5 uM
#* ETO 10 uM

#Time-points:
  
 # * 0h
#* 4h
#* 24h

#Treatment in columns:
  
 # * DMSO: 3
#* ETO 1 uM: 4
#* ETO 2.5 uM: 5
#* ETO 5 uM: 6
#* ETO 10 uM: 7

#Wells in rows:
  
#  * 1: 2
#* 2: 3
#* 3: 4

# get data MDAMB231
data_MDA_0h_1 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT0H-1__2023-10-27T12_38_08-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_0h_2 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT0H-2__2023-10-27T13_34_04-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_2h_1 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT2H-1__2023-10-27T11_30_29-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_2h_2 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT2H-2__2023-10-27T12_09_12-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_4h_1 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT4H-1__2023-10-27T10_50_05-Measurement 1/Evaluation3/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_4h_2 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT4H-2__2023-10-27T11_06_20-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_24h_1 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT24H-1__2023-10-27T10_13_03-Measurement 2/Evaluation1/Objects_Population - Cyto Selected.txt", skip = 9)

data_MDA_24h_2 <- read.delim("/Users/alisc/Downloads/MDA-231/MDA-231/MDAWT24H-2__2023-10-27T10_33_03-Measurement 1/Evaluation1/Objects_Population - Cyto Selected.txt", skip = 9)

# get data Cal51
data_CAL_0h <- read.delim("/Users/alisc/Downloads/cal51.1/cal51.1/T0cal51.imp.parp__2023-03-21T09_12_25-Measurement 1/Evaluation4/Objects_Population - Cyto Selected.txt", skip = 9)

data_CAL_4h <- read.delim("/Users/alisc/Downloads/cal51.1/cal51.1/T4cal51.imp.parp1__2023-03-21T07_39_24-Measurement 1/Evaluation2/Objects_Population - Cyto Selected.txt", skip = 9)

data_CAL_24h <- read.delim("/Users/alisc/Downloads/cal51.1/cal51.1/T24cal51.imp.parp__2023-03-21T10_53_38-Measurement 1/Evaluation3/Objects_Population - Cyto Selected.txt", skip = 9)


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

# lorena wanted me to remove one well: G5 from MDA time 0
data_MDA_all_cond2 <- data_MDA_all_cond %>%
  filter(!(Timepoint == "0h" & Condition == "ETO 2.5uM" &  Well == 3))


# factor
data_MDA_all_cond2$Timepoint <- factor(data_MDA_all_cond2$Timepoint, levels = c("0h", "2h", "4h", "24h"))
data_MDA_all_cond2$Condition <- factor(data_MDA_all_cond2$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_MDA_all_cond2$Well <- factor(data_MDA_all_cond2$Well)
data_MDA_all_cond2$Plate <- factor(data_MDA_all_cond2$Plate)

data_CAL_all_cond$Timepoint <- factor(data_CAL_all_cond$Timepoint, levels = c("0h", "4h", "24h"))
data_CAL_all_cond$Condition <- factor(data_CAL_all_cond$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_CAL_all_cond$Well <- factor(data_CAL_all_cond$Well)

# MDAMB231

# normalise PARP1 cytosolic signal by DMSO

# get DMSO every time point
MDA_0h_DMSO <- data_MDA_all_cond2 %>%
  filter(Timepoint == "0h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

MDA_2h_DMSO <- data_MDA_all_cond2 %>%
  filter(Timepoint == "2h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

MDA_4h_DMSO <- data_MDA_all_cond2 %>%
  filter(Timepoint == "4h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

MDA_24h_DMSO <- data_MDA_all_cond2 %>%
  filter(Timepoint == "24h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

# get median
MDA_0h_DMSO_median <- median(MDA_0h_DMSO, na.rm = T)
MDA_2h_DMSO_median <- median(MDA_2h_DMSO, na.rm = T)
MDA_4h_DMSO_median <- median(MDA_4h_DMSO, na.rm = T)
MDA_24h_DMSOmedian <- median(MDA_24h_DMSO, na.rm = T)

# get Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean normalised by DMSO
data_MDA_all_cond_norm <- data_MDA_all_cond2 %>%
  mutate(PARP1_cyt_norm = case_when(
    Timepoint == "0h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/MDA_0h_DMSO_median,
    Timepoint == "2h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/MDA_2h_DMSO_median,
    Timepoint == "4h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/MDA_4h_DMSO_median,
    Timepoint == "24h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/MDA_24h_DMSOmedian)
  )

# plot PARP1 cytosolic normalised by DMSO
ggplot(data_MDA_all_cond_norm, aes(x = Timepoint, y = log2(PARP1_cyt_norm), fill = Condition)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "Time point", y = "log2(PARP1 cytosol mean intensity \n normalised by DMSO)", 
       fill = "") +
  scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                               "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                               "ETO 10uM" = "#CC6600")) +                       
  theme_classic()

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

PARP1cytnorm <- ggboxplot(
  data_MDA_all_cond_norm, x = "Condition", y = "log2(PARP1_cyt_norm)", fill = "Condition",
  facet.by = "Timepoint", nrow = 1) +
  stat_compare_means(
    comparisons = my_comp, 
    label = "p.format", size = 2
  )

PARP1cytnorm +  labs(y ="log2(PARP1 cytosol mean intensity \n normalised by DMSO)") +
  scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                               "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                               "ETO 10uM" = "#CC6600")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=1))

#remove outlier
data_cond_1_24h <- data_MDA_all_cond_norm[data_MDA_all_cond_norm$Timepoint == "24h",]
data_cond_1 <- data_cond_1_24h %>%
  group_by(Condition) %>%
  mutate(median_PARP1 = median(PARP1_cyt_norm, na.rm = TRUE),
         sd_PARP1 = sd(PARP1_cyt_norm, na.rm = TRUE)) %>%
  filter(PARP1_cyt_norm >= median_PARP1 - (3 * sd_PARP1) &
           PARP1_cyt_norm <= median_PARP1 + (3 * sd_PARP1)) %>%
  ungroup()

pdf("/Users/alisc/Desktop/CRG/final_figures/PARP1_cytosol_MDA.pdf")
ggplot(data_cond_1, aes(Condition, log2(PARP1_cyt_norm), fill= Condition)) + geom_boxplot() +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
   labs(y ="log2(PARP1 cytosol mean intensity \n normalised by DMSO)") +
  scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                               "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                               "ETO 10uM" = "#CC6600")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=1))
dev.off()
table(data_cond_1$Condition)
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


# normalise cytosolic PARP1 by DMSO condition

# get DMSO every time point
CAL_0h_DMSO <- data_CAL_all_cond %>%
  filter(Timepoint == "0h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

CAL_4h_DMSO <- data_CAL_all_cond %>%
  filter(Timepoint == "4h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

CAL_24h_DMSO <- data_CAL_all_cond %>%
  filter(Timepoint == "24h" & Condition == "DMSO") %>%
  pull(Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean)

# get median
CAL_0h_DMSO_median <- median(CAL_0h_DMSO, na.rm = T)
CAL_4h_DMSO_median <- median(CAL_4h_DMSO, na.rm = T)
CAL_24h_DMSO_median <- median(CAL_24h_DMSO, na.rm = T)

# get Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean normalised by DMSO
data_CAL_all_cond_norm <- data_CAL_all_cond %>%
  mutate(PARP1_cyt_norm = case_when(
    Timepoint == "0h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/CAL_0h_DMSO_median,
    Timepoint == "4h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/CAL_4h_DMSO_median,
    Timepoint == "24h" ~ Cyto.Selected...Intensity.Cytoplasm.Alexa.555.Mean/CAL_24h_DMSO_median)
  )


# plot PARP1 cytosol normalised
ggplot(data_CAL_all_cond_norm, aes(x = Timepoint, y = log2(PARP1_cyt_norm), fill = Condition)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "Time point", y = "log2(PARP1 cytosol mean intensity \n normalised by DMSO)", 
       fill = "") +
  scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                               "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                               "ETO 10uM" = "#CC6600")) +                       
  theme_classic()
#ggsave("plots/log2PARP1cyt_CAL_ETO_normDMSO.pdf", device = "pdf", width = 7, height = 4)

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

calparpcytnorm <- ggboxplot(
  data_CAL_all_cond_norm, x = "Condition", y = "log2(PARP1_cyt_norm)", fill = "Condition",
  facet.by = "Timepoint", nrow = 1) +
  stat_compare_means(
    comparisons = my_comp, 
    label = "p.format", size = 2
  )

calparpcytnorm +  labs(y ="log2(PARP1 cytosol mean intensity \n normalised by DMSO)") +
  scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                               "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                               "ETO 10uM" = "#CC6600")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=1))

#remove outlier
data_cond_1_24h <- data_CAL_all_cond_norm[data_CAL_all_cond_norm$Timepoint == "24h",]
data_cond_1 <- data_cond_1_24h %>%
  group_by(Condition) %>%
  mutate(median_PARP1 = median(PARP1_cyt_norm, na.rm = TRUE),
         sd_PARP1 = sd(PARP1_cyt_norm, na.rm = TRUE)) %>%
  filter(PARP1_cyt_norm >= median_PARP1 - (3 * sd_PARP1) &
           PARP1_cyt_norm <= median_PARP1 + (3 * sd_PARP1)) %>%
  ungroup()

pdf("/Users/alisc/Desktop/CRG/final_figures/PARP1_cytosol_CAL.pdf")
ggplot(data_cond_1, aes(Condition, log2(PARP1_cyt_norm), fill= Condition)) + geom_boxplot() +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  labs(y ="log2(PARP1 cytosol mean intensity \n normalised by DMSO)") +
  scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266", 
                               "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000", 
                               "ETO 10uM" = "#CC6600")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=1))
dev.off()
table(data_cond_1$Condition)
