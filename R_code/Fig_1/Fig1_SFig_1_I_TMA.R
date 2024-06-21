library(tidyverse)
library(ggpubr)
library(ggplot2)

# TMA Analysis
#TMA information: https://www.tissuearray.com/tissue-arrays/Breast?product_id=59672 This info is stored in the file TMA_characteristics.
#TMA raw data: 230314_measurements_Eva_Marta_Raul.txt. I have changed the column names so it can be properly read into R into the file: 230314_measurements_Eva_Marta_Raul_NPL.csv. 
#TMA summary data: 230314_TMA_summary.

# get input files
TMA_info <- read_csv("/Users/alisc/Desktop/CRG/final_figures/TMA/TMA_characteristics.csv")
TMA_summary <- read_csv("/Users/alisc/Desktop/CRG/final_figures/TMA/230314_TMA_summary.csv")
TMA_raw <- read_csv("/Users/alisc/Desktop/CRG/final_figures/TMA/230314_measurements_Eva_Marta_Raul_NPL.csv")

# edit positions
colnames(TMA_info)[1] <- "Position"
TMA_summary$`TMA core` <- str_remove(TMA_summary$`TMA core`, "-")
colnames(TMA_summary)[1] <- "Position"
table(table(TMA_summary$Position) == table(TMA_info$Position))
TMA_raw$`TMA core` <- str_remove(TMA_raw$`TMA core`, "-")
colnames(TMA_raw)[5] <- "Position"

# remove position NA
nrow(TMA_raw)
TMA_raw <- TMA_raw %>%
  filter(!is.na(Position))
nrow(TMA_raw) # 684157

# Start with summary data
#Plot the % IMPDH2 nuclear positive cells in:
  # * ER/PR, HER2, TNBC
# add Class, Grade and Tumor stage, Ki67
TMA_info_int <- TMA_info %>%
  dplyr::select(Position, Class, Grade, Stage, Ki67_1, Ki67_2_perc)
TMA_summary_info <- TMA_summary %>%
  left_join(TMA_info_int)

# factor class
TMA_summary_info$Class <- factor(TMA_summary_info$Class, levels = c("ER/PR", "HER2", "TN"))
# remove data with no class
TMA_summary_info_class <- TMA_summary_info %>%
  filter(!is.na(Class))
# remove data with no grade
TMA_summary_info_grade <- TMA_summary_info %>%
  filter(!is.na(Grade))

# remove data with unknown meaning
TMA_summary_info_grade_2 <- TMA_summary_info_grade %>%
  filter(Grade %in% c(1, 2, 3))

# plot % of nuclear positive IMPDH2 cells (Fig_1I)
my_comp <- list(c(1, 2), c(1, 3), c(2, 3))

ggplot(TMA_summary_info_grade_2, aes(x=Grade, y=Percent_positive, fill = Grade)) +
  geom_boxplot(fatten = 2, width = 0.75) +
  labs(x = "Grade", y = "IMPDH2 nuclear positive cells (%)", fill = "") +
  stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
  theme_classic() 

# classify Ki67
TMA_summary_info$Ki67_2_perc <- as.numeric(TMA_summary_info$Ki67_2_perc)

TMA_summary_info_Ki67 <- TMA_summary_info %>%
  filter(Ki67_1 != "*") %>%
  mutate(Ki67_level = case_when(
    Ki67_1 == "-" ~ "Low", 
    Ki67_1 == "+" & Ki67_2_perc < 15 ~ "Low", 
    Ki67_1 == "+" & Ki67_2_perc >= 15 ~ "High"))
TMA_summary_info_Ki67$Ki67_level <- factor(TMA_summary_info_Ki67$Ki67_level, levels = c("Low", "High"))

# From raw data
#calculate from the raw data, the mean/median of nuclear/cytosolic IMPDH2 intensity of all the cells in a specific TMA.
#Variables: Nucleus_DAB_OD_mean and Cytoplasm_DAB_OD_mean
# check type of variable
typeof(TMA_raw$Nucleus_DAB_OD_mean)
typeof(TMA_raw$Cytoplasm_DAB_OD_mean)

# check all TMA
length(sort(unique(TMA_raw$Position)))
sort(unique(TMA_raw$Position))

# get the mean/median values of nuclear IMPDH2 intensity per position (TMA)
TMA_raw_summary <- TMA_raw %>%
  group_by(Position) %>%
  summarise(n = n(),
            mean_IMPDH2_nucleus = mean(Nucleus_DAB_OD_mean),
            median_IMPDH2_nucleus = median(Nucleus_DAB_OD_mean),
            stdev_IMPDH2_nucleus = sd(Nucleus_DAB_OD_mean),
            se_IMPDH2_nucleus = stdev_IMPDH2_nucleus/sqrt(n),
            mean_IMPDH2_cyt = mean(Cytoplasm_DAB_OD_mean),
            median_IMPDH2_cyt = median(Cytoplasm_DAB_OD_mean),
            stdev_IMPDH2_cyt = sd(Cytoplasm_DAB_OD_mean),
            se_IMPDH2_cyt = stdev_IMPDH2_cyt/sqrt(n))

# add Class, Grade and Tumor stage, Ki67
TMA_raw_summary_info <- TMA_raw_summary %>%
  left_join(TMA_info_int)

# factor class
TMA_raw_summary_info$Class <- factor(TMA_raw_summary_info$Class, levels = c("ER/PR", "HER2", "TN"))

# remove data with no class
TMA_raw_summary_infoclass <- TMA_raw_summary_info %>%
  filter(!is.na(Class))

# add new variable
TMA_raw_summary_infoclass_2 <- TMA_raw_summary_infoclass %>%
  mutate(TNBC = ifelse(Class == "TN", "TN", "non-TN"))
table(TMA_raw_summary_infoclass_2$TNBC)
TMA_raw_summary_infoclass_2$TNBC <- factor(TMA_raw_summary_infoclass_2$TNBC, levels = c("non-TN", "TN"))

#Fig 1 I part 2
# remove data with no grade
TMA_raw_summary_info_grade <- TMA_raw_summary_info %>%
  filter(!is.na(Grade))

table(TMA_raw_summary_info_grade$Grade)

# remove data with unknown meaning
TMA_raw_summary_info_grade_2 <- TMA_raw_summary_info_grade %>%
  filter(Grade %in% c(1, 2, 3))

table(TMA_raw_summary_info_grade_2$Grade)

# plot IMPDH2 nuclear intensity: mean and median
my_comp <- list(c(1, 2), c(1, 3), c(2, 3))

ggplot(TMA_raw_summary_info_grade_2, aes(x=Grade, y=mean_IMPDH2_nucleus, fill = Grade)) +
  geom_boxplot(fatten = 2, width = 0.75) +
  labs(x = "Grade", y = "IMPDH2 mean nuclear signal", fill = "") +
  stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
  theme_classic() 

#Supplementary Fig 1I
# plot IMPDH2 nuclear intensity: mean and median
my_comp <- list(c("non-TN", "TN"))

ggplot(TMA_raw_summary_infoclass_2, aes(x=TNBC, y=mean_IMPDH2_nucleus, fill = TNBC)) +
  geom_boxplot(fatten = 2, width = 0.75) +
  labs(x = "", y = "IMPDH2 mean nuclear signal", fill = "") +
  stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
  scale_fill_manual(name = "", values = c("non-TN" = "#CAD2C5", "TN" = "#354f52"))+
  theme_classic() 

#remove negative values
TMA_raw_summary_infoclass_2_no_neg <- TMA_raw_summary_infoclass_2[TMA_raw_summary_infoclass_2$mean_IMPDH2_nucleus > 0,]

pdf("/Users/alisc/Desktop/CRG/final_figures/TMA/IMPDH2_signal_no_neg.pdf")
ggplot(TMA_raw_summary_infoclass_2_no_neg, aes(x=TNBC, y=mean_IMPDH2_nucleus, fill = TNBC)) +
  geom_boxplot(fatten = 2, width = 0.75) +
  labs(x = "", y = "IMPDH2 mean nuclear signal", fill = "") +
  stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
  scale_fill_manual(name = "", values = c("non-TN" = "#CAD2C5", "TN" = "#354f52"))+
  theme_classic() 
dev.off()
table(TMA_raw_summary_infoclass_2_no_neg$TNBC)
#non-TN     TN 
#66     25 
