library(ggplot2)
library(ggpubr)
library(dplyr)

#load data
#row 2,3,4 are WT
#row 5,6,7 are KO

data <- read.delim("/Users/alisc/Desktop/CRG/Operetta/Guanosine_Titration_27_02_24/Objects_Population - full nuclei.txt", sep="\t")

#Add concentrations for the rows
data$Guanosine <- NA  
data[data$Column == 2, "Guanosine"] <- 0
data[data$Column == 3, "Guanosine"] <- 1000
data[data$Column == 4, "Guanosine"] <- 800
data[data$Column == 5, "Guanosine"] <- 400
data[data$Column == 6, "Guanosine"] <- 200
data[data$Column == 7, "Guanosine"] <- 100
data[data$Column == 8, "Guanosine"] <- 50
data[data$Column == 9, "Guanosine"] <- 25
data[data$Column == 10, "Guanosine"] <- 12.5
data[data$Column == 11, "Guanosine"] <- 6.25

WT_row <- c(2,3,4)
KO_row <- c(5,6,7)
data_WT <- data[data$Row %in% WT_row,]
data_KO <- data[data$Row %in% KO_row,]


ggplot(data_WT, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots, fill=factor(Guanosine))) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of γH2AX foci per nucleus") +
  theme_minimal() +
  ggtitle("γH2AX foci per cell in WT") + theme(legend.position = "none") + scale_fill_grey(start = 1, end = 0.5)+ coord_cartesian(ylim = c(NA, 260))

ggplot(data_KO, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots, fill=factor(Guanosine))) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of γH2AX foci per nucleus") +
  theme_minimal() +
  ggtitle("γH2AX foci per cell in KO") + theme(legend.position = "none")   + scale_fill_grey(start = 1, end = 0.5) + coord_cartesian(ylim = c(NA, 260))


#plot 53BP1
ggplot(data_WT, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots..2., fill=factor(Guanosine))) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of 53BP1 foci per nucleus") +
  theme_minimal() + coord_cartesian(ylim = c(NA, 200))+
  ggtitle("53BP1 foci per cell in WT") + theme(legend.position = "none")+ scale_fill_grey(start = 1, end = 0.5)

ggplot(data_KO, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots..2., fill=factor(Guanosine))) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of 53BP1 foci per nucleus") +
  theme_minimal() + 
  coord_cartesian(ylim = c(NA, 200)) +
  ggtitle("53BP1 foci per cell in KO") + theme(legend.position = "none")+ scale_fill_grey(start = 1, end = 0.5)

#cell count
data$Cell <- ifelse(data$Row %in% WT_row,"WT", "KO" )
summary_count <- data %>%
  group_by(Guanosine,Cell ) %>%
  summarise(count = n())

ggplot(summary_count, aes(factor(Guanosine), count, fill=Cell)) +
  geom_bar(stat="identity", position = "dodge", color="black", width = 0.5) + 
  theme_minimal() +  
  labs(y = "cell count", title = "Cell count") + 
scale_fill_manual(values = c("KO" = "black", "WT" = "grey")) + 
  labs(x = "Guanosine Concentration", y = "Number of cells")

#scale to percentage WT basal
wt_basal <- 10683
summary_count$percentage <- summary_count$count / wt_basal * 100

ggplot(summary_count, aes(factor(Guanosine), percentage, fill=Cell)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width = 0.5) + 
  theme_minimal() +  labs(y = "cell count", title = "Cell count")  + 
  scale_fill_manual(values = c("KO" = "black", "WT" = "grey")) + 
  labs(x = "Guanosine Concentration", y = "percentage of WT basal")

#yH2AX and 53BP1 together
data$Type <- ifelse(data$Row %in% WT_row, "WT", "KO")

ggplot(data, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots..2., fill=Type)) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of 53BP1 foci per nucleus") +
  theme_minimal() +
  ggtitle("53BP1 foci per cell") + scale_fill_grey(start = 1, end = 0.5) 

ggplot(data, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots, fill=Type)) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of yH2AX foci per nucleus") +
  theme_minimal() +
  ggtitle("yH2AX foci per cell") + scale_fill_grey(start = 1, end = 0.5) 

#remove outliers with 5sd
# Add SD and mean for 53BP1
data_53BP1 <- data %>%
  group_by(Guanosine) %>%
  mutate(
    median_53BP1 = median(full.nuclei...Number.of.Spots..2.),
    sd_53BP1 = sd(full.nuclei...Number.of.Spots..2.)
  ) %>%
  ungroup()

filtered_data_53P1 <- data_53BP1 %>%
  filter(
    full.nuclei...Number.of.Spots..2. >= (median_53BP1 - 5 * sd_53BP1) &
      full.nuclei...Number.of.Spots..2. <= (median_53BP1 + 5 * sd_53BP1)
  )

# Add SD and mean for gH2X
data_yH2AX <- data %>%
  group_by(Guanosine) %>%
  mutate(
    median_gH2 = median(full.nuclei...Number.of.Spots),
    sd_gH2 = sd(full.nuclei...Number.of.Spots..2.)
  ) %>%
  ungroup()

filtered_data_yH2AX <- data_yH2AX %>%
  filter(
    full.nuclei...Number.of.Spots >= (median_gH2 - 5 * sd_gH2) &
      full.nuclei...Number.of.Spots <= (median_gH2 + 5 * sd_gH2)
  )

ggplot(filtered_data_yH2AX, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots..2., fill=Type)) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of 53BP1 foci per nucleus") +
  theme_minimal() +
  ggtitle("53BP1 foci per cell") + scale_fill_grey(start = 1, end = 0.5) 

ggplot(filtered_data_53P1, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots, fill=Type)) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of yH2AX foci per nucleus") +
  theme_minimal() +
  ggtitle("yH2AX foci per cell") + scale_fill_grey(start = 1, end = 0.5) 

#remove outliers with +/- 3sd
filtered_data_53P1_2 <- data_53BP1 %>%
  filter(
    full.nuclei...Number.of.Spots..2. >= (median_53BP1 - 3 * sd_53BP1) &
      full.nuclei...Number.of.Spots..2. <= (median_53BP1 + 3 * sd_53BP1)
  )

filtered_data_yH2AX_2 <- data_yH2AX %>%
  filter(
    full.nuclei...Number.of.Spots >= (median_gH2 - 3 * sd_gH2) &
      full.nuclei...Number.of.Spots <= (median_gH2 + 3 * sd_gH2)
  )
#plot
ggplot(filtered_data_yH2AX_2, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots..2., fill=Type)) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of 53BP1 foci per nucleus") +
  theme_minimal() +
  ggtitle("53BP1 foci per cell") + scale_fill_grey(start = 1, end = 0.5) 

ggplot(filtered_data_53P1_2, aes(x=factor(Guanosine), y = full.nuclei...Number.of.Spots, fill=Type)) +
  geom_boxplot() + 
  labs(x = "Guanosine Concentration", y = "Number of yH2AX foci per nucleus") +
  theme_minimal() +
  ggtitle("yH2AX foci per cell") + scale_fill_grey(start = 1, end = 0.5) 
