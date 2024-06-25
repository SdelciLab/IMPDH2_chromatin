library(tidyverse)
library(ggpubr)
library(ggplot2)

#load data
#row 2: DMSO
#row 3: 10uM ATMi
#row 4: 10uM ATRi
#column 2,3,4: 10uM Etoposide
#column 6,7,8: DMSO

data <- read.delim("/Users/alisc/Desktop/CRG/Operetta/ATM_ATR_inhibition/ATM_ATR_plate2/Objects_Population - full nuclei.txt", skip = 9)

#label DMSO and etoposide
Etoposide <- c(2,3,4)
data_cond <- data
data_cond$Treatment <- ifelse(data$Column %in% Etoposide,"Etoposide", "DMSO" )

#label drug 
data_cond$Drug <- NA
data_cond[data_cond$Row == 2, "Drug"] <- "DMSO"
data_cond[data_cond$Row == 3, "Drug"] <- "ATMi"
data_cond[data_cond$Row == 4, "Drug"] <- "ATRi"
data_cond$Identifier <- paste(data_cond$Treatment, data_cond$Drug)

#Filtering
ggplot(data_cond, aes(full.nuclei....Roundness)) + geom_density()
median(data_cond$full.nuclei....Roundness)
# most cells have a roundness of 0.9123, filtering for min 0.8 roundness 
filtered_data <- data_cond[data_cond$full.nuclei....Roundness > 0.8, ]
ggplot(filtered_data, aes(full.nuclei....Roundness)) + geom_density()

ggplot(data_cond, aes(full.nuclei....Area..µm..)) + geom_density()
filtered_data2 <- filtered_data[filtered_data$full.nuclei....Area..µm.. < 400, ]

#plot yH2AX
filtered_data2$Identifier <- factor(filtered_data2$Identifier, 
                                    levels = c("DMSO DMSO", "DMSO ATMi","DMSO ATRi","Etoposide DMSO","Etoposide ATMi", "Etoposide ATRi"))
my_comparisons <- list(c("DMSO DMSO", "DMSO ATMi"), c("DMSO DMSO", "DMSO ATRi"), c("Etoposide DMSO", "Etoposide ATMi"),  c("Etoposide DMSO", "Etoposide ATRi"))

ggplot(filtered_data2, aes(x = Identifier, y =full.nuclei...Number.of.Spots , fill = Treatment)) + 
  geom_boxplot(width = 0.75, fatten = 2)  +   
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.signif") +
  labs(x = "", y = "yH2AX foci per nucleus") + 
  theme_classic()+ scale_fill_brewer(palette="Accent") +  theme(legend.position = "none")

ggplot(filtered_data2, aes(x = Drug, y =full.nuclei...Number.of.Spots , fill = Treatment)) +
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "yH2AX foci per nucleus") +
  theme_classic() +  scale_fill_brewer(palette="Accent") + 
  stat_compare_means(method = "wilcox", aes(group = Treatment), label = "p.signif")  

#plot IMPDH2
ggplot(filtered_data2, aes(x = Identifier, y = full.nuclei...Intensity.Nucleus.Alexa.488.Mean, fill= Treatment)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "IMPDH2 intensity nucleus") +  
  theme_classic() +  scale_fill_brewer(palette="Accent")+  theme(legend.position = "none") +   
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.signif")

ggplot(filtered_data2, aes(x = Drug, y = full.nuclei...Intensity.Nucleus.Alexa.488.Mean, fill= Treatment)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "IMPDH2 intensity nucleus")  +
  stat_compare_means(method = "wilcox", aes(group = Treatment), label = "p.signif") +
  theme_classic()+  scale_fill_brewer(palette="Accent")

#plot yH2AX normalised by DMSO
filtered_data2 <- filtered_data2 %>%
  group_by(Drug, Treatment) %>%
  mutate(Median_Value = median(full.nuclei...Number.of.Spots, na.rm = TRUE))

medians_DMSO <- filtered_data2 %>% 
  filter(Treatment == "DMSO") %>%
  group_by(Drug) %>%
  summarise(Median_DMSO = median(full.nuclei...Number.of.Spots, na.rm = TRUE))

# Step 2: Merge this median back to the original dataframe and normalize the values
df_normalized <- filtered_data2 %>%
  left_join(medians_DMSO, by = "Drug") %>%
  mutate(Normalized_Value = full.nuclei...Number.of.Spots - Median_DMSO)

ggplot(df_normalized, aes(x = Drug, y =Normalized_Value , fill = Treatment)) +
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "yH2AX foci normalised by DMSO") +
  theme_classic() +  scale_fill_brewer(palette="Accent") + 
  stat_compare_means(method = "wilcox", aes(group = Treatment), label = "p.signif") 

#plot IMPDH2 normalised by DMSO
filtered_data2 <- filtered_data2 %>%
  group_by(Drug, Treatment) %>%
  mutate(Median_Value = median(full.nuclei...Intensity.Nucleus.Alexa.488.Mean, na.rm = TRUE))

medians_DMSO <- filtered_data2 %>% 
  filter(Treatment == "DMSO") %>%
  group_by(Drug) %>%
  summarise(Median_DMSO = median(full.nuclei...Intensity.Nucleus.Alexa.488.Mean, na.rm = TRUE))

# Step 2: Merge this median back to the original dataframe and normalize the values
df_normalized <- filtered_data2 %>%
  left_join(medians_DMSO, by = "Drug") %>%
  mutate(Normalized_Value = full.nuclei...Intensity.Nucleus.Alexa.488.Mean - Median_DMSO)

ggplot(df_normalized, aes(x = Drug, y =Normalized_Value , fill = Treatment)) +
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "IMPDH2 intensity normalised by DMSO") +
  theme_classic() +  scale_fill_brewer(palette="Accent") + 
  stat_compare_means(method = "wilcox", aes(group = Treatment), label = "p.signif")  

