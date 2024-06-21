library(tidyverse)
library(ggpubr)
library(pracma)

#Row of interest
#Row 2: MCF7
#Row 3: MCFZ hormone insensitive (HI)
#Row 4: MDAMB231 (not include)

#load data
data <- read.delim("/Users/alisc/Desktop/CRG/Operetta/HIplate_Marta/DAPI/Objects_Population - Nuclei Selected 5.txt", skip = 9)
data_cond <- data %>%
  mutate(Cell = case_when(
    Row == 2 ~ "MCF7",
    Row == 3 ~ "HI",
    Row == 4 ~ "MDAMB231")
  )
data_cond$Cell <- factor(data_cond$Cell, levels = c("MCF7", "HI","MDAMB231"))
data_filtered <- data_cond %>%
  filter(Cell %in% c("HI", "MCF7"))
my_comp <- list(c("HI", "MCF7"))

#add corrected DAPI signal by substracting RING region intensity
data_filtered$DAPIcorrected <- data_filtered$Nuclei.Selected...Intensity.Nucleus.DAPI.Mean - data_filtered$Nuclei.Selected...Intensity.Ring.Region.DAPI.Mean

#calculate inttegrated signal intensity
data_filtered$ratio <- data_filtered$DAPIcorrected*data_filtered$Nuclei.Selected...Nucleus.Area..µm..
data_filtered <- data_filtered[data_filtered$ratio<1e+07,] #remove outliers

#Normalise G1 peaks
a1 <- ggplot(data_filtered, aes(ratio, color=Cell)) +geom_density()
density_data2 <- ggplot_build(a1)$data[[1]]
group_levels2 <- levels(factor(data_filtered$Cell))
density_data2$Cell <- group_levels2[density_data2$group]

#check for peak value per condtion
density2 <- density_data2 %>%
  group_by(Cell) %>%
  mutate(MaxY = max(y)) %>% 
  filter(y == MaxY) %>% 
  summarise(x = first(x), MaxY = first(MaxY), .groups = 'drop')

# WT peak x-coordinate
HI_peak_x <- density2$x[density2$Cell == "HI"]

# Calculate shifts needed for each condition
density2$shift2 = HI_peak_x - density2$x

# Join data with shifts
normalised_data_G1 <- data_filtered %>%
  left_join(density2, by = "Cell") %>%
  mutate(adjusted_value_G1 = ratio + shift2) 

#check cell cycle plot after normalising G1 peaks
ggplot(normalised_data_G1, aes(adjusted_value_G1, color=Cell)) +geom_density() 

#subset data for proliferating cells only (G1/G0 vs rest)
#Extract the values of the peaks and valleys
data_HI <- normalised_data_G1 %>%
  dplyr::filter(Cell %in% c("HI"))

dens <- density(data_HI$adjusted_value_G1)

peaks <- findpeaks(dens$y)
peaks <- as.data.frame(peaks)
peak_x_values <- dens$x[peaks$V2]

troughs <- findpeaks(-dens$y)
troughs <- as.data.frame(troughs)
troughs_x_values <- dens$x[troughs$V2]
#come up with gates
gate1 <- peak_x_values[4] -500000
gate2 <- peak_x_values[4] + 500000
gate3 <- troughs_x_values[4] -500000
gate4 <- troughs_x_values[4] +500000

ggplot(normalised_data_G1, aes(adjusted_value_G1, color=Cell)) + geom_density() +  # Density curve
  # geom_vline(xintercept = peak_x_values[4], col = "red", linetype = "dashed") +
  geom_vline(xintercept = troughs_x_values[4], col = "red", linetype = "dashed") + 
  labs(title = "Density with Peaks Marked", x = "DAPI integrated intensity", y = "") +
  theme_bw() + scale_color_brewer(palette="Accent")

#subset and balance data
subset_data <- normalised_data_G1 %>%
  filter(adjusted_value_G1 > troughs_x_values[4])

table(subset_data$Cell)
data_HI <- subset_data %>%
  filter(Cell == "HI")
seed <- set.seed(200)
# Filter data for MCF7 and then sample 363  data points
data_MCF7 <- subset_data %>%
  filter(Cell == "MCF7") %>%
  sample_n(size = 363 )

# Combine the two balanced datasets
balanced_data <- bind_rows(data_HI, data_MCF7)
table(balanced_data$Cell)

#plot yH2AX foci (Nuclei.Selected...Number.of.Spots)
my_comp <- list(c("HI", "MCF7"))

ggplot(balanced_data, aes(x = Cell, y = Nuclei.Selected...Number.of.Spots, color = Cell)) + geom_jitter(alpha=0.4,  width = 0.2) +
  geom_boxplot(color= "black", fill= NA, width = 0.5) +
  labs(x = "", y = "Number of yH2AX spots per nucleus", fill = "") +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() + scale_color_brewer(palette="Accent") +
  theme( legend.position = "none")
