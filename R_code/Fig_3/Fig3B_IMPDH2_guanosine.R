library(tidyverse)
library(ggpubr)
library(ggplot2)

#Compare gH2AX foci and IMPDH2 signal in nucleus
#Row 2: DMSO control
#Row 3: NAD+ (100uM per well)
#Row 4: Guanosine (50uM per well)
#Column 2,3,4: DMSO control
#Column 6,7,8: Etoposide treatment 3h (10uM per well)

data <- read.delim("/Users/alisc/Desktop/CRG/Operetta/NADGuanosineplate_1_and_2/Objects_Population - full nuclei_Guanosine_IMPDH2.txt", skip = 9)

#label DMSO and etoposide
DMSO <- c(2,3,4)
data_cond <- data
data_cond$Condition <- ifelse(data$Column %in% DMSO,"DMSO", "Etoposide")

#label treatment 
data_cond$treatment <- NA
data_cond[data_cond$Row == 2, "treatment"] <- "DMSO"
data_cond[data_cond$Row == 3, "treatment"] <- "NAD"
data_cond[data_cond$Row == 4, "treatment"] <- "Guanosine"

data_cond$Identifier <- paste(data_cond$Condition, data_cond$treatment)
data_cond$Identifier <- factor(data_cond$Identifier, 
                               levels = c("DMSO DMSO","DMSO Guanosine", "DMSO NAD","Etoposide DMSO","Etoposide Guanosine", "Etoposide NAD"))
data_cond$treatment <- factor(data_cond$treatment, 
                              levels = c("DMSO","Guanosine", "NAD"))
#subset data without NAD+
noNAD <-data_cond[data_cond$treatment != "NAD",]

#plt gH2AX signal
ggplot(noNAD, aes(x = treatment, y =full.nuclei...Number.of.Spots, fill = Condition)) + 
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "gH2AX foci per nucleus") + theme(legend.position = "none")+
  theme_classic() +
  scale_fill_brewer(palette="Accent")

ggplot(noNAD, aes(x = Identifier, y =full.nuclei...Number.of.Spots , fill = Condition)) + 
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "gH2AX foci per nucleus") +
  theme_classic()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme(legend.position = "none")+
  scale_fill_brewer(palette="Accent")

#plot log IMPDH2 signal
my_comparisons <- list(c("DMSO DMSO", "DMSO Guanosine"), c("Etoposide DMSO", "Etoposide Guanosine"))

ggplot(noNAD, aes(x = treatment, y =log1p(full.nuclei...Intensity.Nucleus.Alexa.488.Mean), fill = Condition)) +
  geom_boxplot(width = 0.75, fatten = 2) + stat_compare_means(method = "wilcox", aes(group = Condition), label = "p.signif") +
  labs(x = "", y = "IMPDH2 mean intensity nucleus") +
  theme_classic() +
  scale_fill_brewer(palette="Accent")  + theme(legend.position = "none")
ggplot(noNAD, aes(x = Identifier, y = log1p(full.nuclei...Intensity.Nucleus.Alexa.488.Mean), fill= Condition)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "IMPDH2 mean intensity nucleus") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme(legend.position = "none")  + stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.signif")+
  scale_fill_brewer(palette="Accent")

#remove outlier
noNAD_outlier <- noNAD %>%
  
  group_by(Identifier) %>%
  mutate(median_1 = median(full.nuclei...Intensity.Nucleus.Alexa.488.Mean),
         sd_1 = sd(full.nuclei...Intensity.Nucleus.Alexa.488.Mean)) %>%
  
  filter(full.nuclei...Intensity.Nucleus.Alexa.488.Mean >= median_1 - (3 * sd_1) &
           full.nuclei...Intensity.Nucleus.Alexa.488.Mean <= median_1 + (3 * sd_1)) %>%
  ungroup() 
pdf("/Users/alisc/Desktop/CRG/final_figures/NAD_Guanosine1/Guanosine_noOutlier.pdf")
ggplot(noNAD_outlier, aes(x = Identifier, y = log1p(full.nuclei...Intensity.Nucleus.Alexa.488.Mean), fill= Condition)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "IMPDH2 mean intensity nucleus") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme(legend.position = "none")  + stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  scale_fill_brewer(palette="Accent")
dev.off()