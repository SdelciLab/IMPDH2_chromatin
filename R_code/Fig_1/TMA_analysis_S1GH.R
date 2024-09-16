# libraries
library(tidyverse)
library(ggpubr)

# get input files
TMA_info <- read_csv("input data/TMA_characteristics.csv")
TMA_summary <- read_csv("input data/230314_TMA_summary.csv")
TMA_raw <- read_csv("input data/230314_measurements_Eva_Marta_Raul_NPL.csv")

# explore TMA info
head(TMA_info)
str(TMA_info)
colnames(TMA_info)
table(TMA_info$Class)
table(is.na(TMA_info$Class))
table(TMA_info$Stage)
table(TMA_info$Grade)
colnames(TMA_info)[1] <- "Position"

# check all positions
table(TMA_info$Position)

# explore summary
head(TMA_summary)
str(TMA_summary)
colnames(TMA_summary)

# edit position 
TMA_summary$`TMA core` <- str_remove(TMA_summary$`TMA core`, "-")
colnames(TMA_summary)[1] <- "Position"
table(table(TMA_summary$Position) == table(TMA_info$Position))

# explore TMA raw data
head(TMA_raw)
str(TMA_raw)
colnames(TMA_raw)

# edit position 
TMA_raw$`TMA core` <- str_remove(TMA_raw$`TMA core`, "-")
colnames(TMA_raw)[5] <- "Position"

# check
table(table(TMA_summary$Position) == table(TMA_info$Position)) # TRUE
table(table(TMA_summary$Position) == table(TMA_raw$Position)) # False 
table(sort(unique(TMA_summary$Position)) == sort(unique(TMA_raw$Position))) # TRUE

# check na
table(TMA_raw$Position)
sum(table(TMA_raw$Position)) # 684157
table(is.na(TMA_raw$Position)) # True=731 

# remove position NA
nrow(TMA_raw)
TMA_raw <- TMA_raw %>%
    filter(!is.na(Position))
nrow(TMA_raw) # 684157

## SUMMARY DATA

# add Class, Grade and Tumor stage, Ki67
TMA_info_int <- TMA_info %>%
    dplyr::select(Position, Class, Grade, Stage, Ki67_1, Ki67_2_perc)

TMA_summary_info <- TMA_summary %>%
    left_join(TMA_info_int)

# remove data with no grade
TMA_summary_info_grade <- TMA_summary_info %>%
    filter(!is.na(Grade))

table(TMA_summary_info_grade$Grade)

# remove data with unknown meaning
TMA_summary_info_grade_2 <- TMA_summary_info_grade %>%
    filter(Grade %in% c(1, 2, 3))

table(TMA_summary_info_grade_2$Grade)

# plot % of nuclear positive IMPDH2 cells
my_comp <- list(c(1, 2), c(1, 3), c(2, 3))

ggplot(TMA_summary_info_grade_2, aes(x=Grade, y=Percent_positive, fill = Grade)) +
    geom_boxplot(fatten = 2, width = 0.75) +
    labs(x = "Grade", y = "IMPDH2 nuclear positive cells (%)", fill = "") +
    stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
    theme_classic() 
ggsave("plots_summary/perc_IMPDH2_nuc_grade.pdf", device = "pdf", width = 5, height = 4)

## RAW DATA

# check type of variable
typeof(TMA_raw$Nucleus_DAB_OD_mean)
typeof(TMA_raw$Cytoplasm_DAB_OD_mean)

# check all TMA
length(sort(unique(TMA_raw$Position)))
sort(unique(TMA_raw$Position))
table(sort(unique(TMA_raw$Position)) == sort(unique(TMA_info$Position))) # T

# check na
nrow(TMA_raw) # 684157
table(is.na(TMA_raw$Cytoplasm_DAB_OD_mean)) # none
table(is.na(TMA_raw$Nucleus_DAB_OD_mean)) # none

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
head(TMA_raw_summary)
nrow(TMA_raw_summary)

# check all TMA
length(sort(unique(TMA_raw_summary$Position)))
sort(unique(TMA_raw_summary$Position))

# add Class, Grade and Tumor stage, Ki67
TMA_raw_summary_info <- TMA_raw_summary %>%
    left_join(TMA_info_int)

# factor class
TMA_raw_summary_info$Class <- factor(TMA_raw_summary_info$Class, levels = c("ER/PR", "HER2", "TN"))

# remove data ith no class
TMA_raw_summary_infoclass <- TMA_raw_summary_info %>%
    filter(!is.na(Class))

table(TMA_raw_summary_infoclass$Class)

# remove data with no grade
TMA_raw_summary_info_grade <- TMA_raw_summary_info %>%s
    filter(!is.na(Grade))

table(TMA_raw_summary_info_grade$Grade)

# remove data with unknown meaning
s <- TMA_raw_summary_info_grade %>%
    filter(Grade %in% c(1, 2, 3))

table(TMA_raw_summary_info_grade_2$Grade)

# plot IMPDH2 nuclear intensity: mean and median
my_comp <- list(c(1, 2), c(1, 3), c(2, 3))

ggplot(TMA_raw_summary_info_grade_2, aes(x=Grade, y=mean_IMPDH2_nucleus, fill = Grade)) +
    geom_boxplot(fatten = 2, width = 0.75) +
    labs(x = "Grade", y = "IMPDH2 mean nuclear signal", fill = "") +
    stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
    theme_classic() 
ggsave("plots_ind/meannuc_IMPDH2_grade.pdf", device = "pdf", width = 5, height = 4)