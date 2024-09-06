library(tidyverse)
library(ggpubr)

#Compare mono/poly ADP ribose signal in the nucleus
#row 2 wt-empty
#row 3 ko-empty
#row 4 ko-wt
#row 5 ko-nls
#row 6 wt
#row 7 ko

#load the data
data <- read.delim("/Users/alisc/Desktop/CRG/Operetta/PARplate/Objects_Population - Nuclei 2.txt", skip = 9)

data_cond <- data %>%
  mutate(Cell = case_when(
    Row == 2 ~ "WT_empty",
    Row == 3 ~ "KO_empty",
    Row == 4 ~ "KO-WT",
    Row == 5 ~ "KO-NLS",
    Row == 6 ~ "WT",
    Row == 7 ~ "KO",
  )
  )

#filtering
filtered_data_new <- data_cond[data_cond$Nuclei...Nucleus.Area< 250, ]
ordered <- c("WT", "KO","WT_empty","KO_empty","KO-WT","KO-NLS")
filtered_data_new$Cell <- factor(filtered_data_new$Cell, levels = ordered)
table(filtered_data_new$Cell)

#Plot mean mono/poly ADP ribose signal nucleus
my_comp <- list(c("WT", "KO"), c("WT_empty", "KO_empty"), c("KO-WT", "KO-NLS"))
ggplot(filtered_data_new, aes(x = Cell, y = Nuclei...Intensity.Nucleus.Alexa.488.Mean, fill = Cell)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "", y = "Intensity Mono/Poly ADP-ribose", fill = "") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() +  scale_fill_brewer(palette="Accent") +
  theme( legend.position = "none")

#seperate WT KO
WT_KO <-  filtered_data_new[filtered_data_new$Row %in% c(6,7),]
recon <- filtered_data_new[filtered_data_new$Row %in% c(2,3,4,5),]
my_comp <- list(c("WT", "KO"))
ggplot(WT_KO, aes(x = Cell, y = Nuclei...Intensity.Nucleus.Alexa.488.Mean, fill = Cell)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "", y = "Intensity Mono/Poly ADP-ribose", fill = "") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() +  scale_fill_brewer(palette="Accent") +
  theme( legend.position = "none")

my_comp <- list(c("WT_empty", "KO_empty"), c("KO-WT", "KO-NLS"))
ggplot(recon, aes(x = Cell, y = Nuclei...Intensity.Nucleus.Alexa.488.Mean, fill = Cell)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "", y = "Intensity Mono/Poly ADP-ribose", fill = "") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() +  scale_fill_brewer(palette="Accent") +
  theme( legend.position = "none")
