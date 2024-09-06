# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

#Conditions:
#* WT-control
#* KO-control
#* KO-wt-IMPDH2
#* KO-NES-IMPDH2
#* KO-3xNLS-IMPDH2

#Plate organisation:
#Cells columns:
#* WT-control: 2
#* KO-control: 3
#* KO-wt-IMPDH2: 4 
#* KO-NES-IMPDH2: 5
#* KO-3xNLS-IMPDH2: 6

#Technical replicates in rows: 
#* 1: 5 (E)
#* 2: 6 (F)
#* 3: 7 (G)
#Biological replicates in different plates

```{r include=FALSE}
# get data
H2AX_quant_rep1 <- read.delim("/Users/alisc/Downloads/Objects_Population - full nuclei_recons_130923.txt", skip = 9)
head(H2AX_quant_rep1)

H2AX_quant_rep2 <- read.delim("/Users/alisc/Downloads/Objects_Population - full nuclei_recons210923_replica2.txt", skip = 9)
head(H2AX_quant_rep2)

H2AX_quant_rep3 <- read.delim("/Users/alisc/Downloads/Objects_Population - full nucleireplica3.txt", skip = 9)
head(H2AX_quant_rep3)

# check
colnames(H2AX_quant_rep1)
table(colnames(H2AX_quant_rep1) == colnames(H2AX_quant_rep2))
table(colnames(H2AX_quant_rep1) == colnames(H2AX_quant_rep3))

# plate distribution
table(H2AX_quant_rep1$Row, H2AX_quant_rep1$Column)
table(H2AX_quant_rep2$Row, H2AX_quant_rep2$Column)
table(H2AX_quant_rep3$Row, H2AX_quant_rep3$Column)

# add replicate
H2AX_quant_rep1 <- H2AX_quant_rep1 %>%
  mutate(Replicate = rep(1, nrow(H2AX_quant_rep1)))
H2AX_quant_rep2 <- H2AX_quant_rep2 %>%
  mutate(Replicate = rep(2, nrow(H2AX_quant_rep2)))
H2AX_quant_rep3 <- H2AX_quant_rep3 %>%
  mutate(Replicate = rep(3, nrow(H2AX_quant_rep3)))

# put together
H2AX_quant_allreps <- H2AX_quant_rep1 %>%
  rbind(H2AX_quant_rep2) %>%
  rbind(H2AX_quant_rep3)

# add condition
H2AX_quant_allreps_cond <- H2AX_quant_allreps %>%
  mutate(Condition = case_when(
    Column == 2 ~ "WT-control",
    Column == 3 ~ "KO-control",
    Column == 4 ~ "KO-wt-IMPDH2",
    Column == 5 ~ "KO-NES-IMPDH2",
    Column == 6 ~ "KO-3xNLS-IMPDH2"),
    Well = case_when(
      (Row == 5) ~ 1,
      (Row == 6) ~ 2,
      (Row == 7) ~ 3)
  )

# filter interesting data
H2AX_quant_allreps_cond_filt <- H2AX_quant_allreps_cond %>%
  filter(Condition %in% c("WT-control", "KO-control", "KO-wt-IMPDH2", "KO-NES-IMPDH2", "KO-3xNLS-IMPDH2")) %>%
  filter(Replicate %in% c(1,2,3)) %>%
  filter(Well %in% c(1,2,3))

# verify
table(H2AX_quant_allreps_cond_filt$Condition, H2AX_quant_allreps_cond_filt$Column)
table(H2AX_quant_allreps_cond_filt$Well, H2AX_quant_allreps_cond_filt$Row)

table(H2AX_quant_allreps_cond_filt$Replicate, H2AX_quant_allreps_cond_filt$Row)
table(H2AX_quant_allreps_cond_filt$Replicate, H2AX_quant_allreps_cond_filt$Column)

# factors
H2AX_quant_allreps_cond_filt$Condition <- factor(H2AX_quant_allreps_cond_filt$Condition, levels = 
                                                   c("WT-control", "KO-control", "KO-wt-IMPDH2", "KO-3xNLS-IMPDH2",
                                                     "KO-NES-IMPDH2"))
H2AX_quant_allreps_cond_filt$Replicate <- factor(H2AX_quant_allreps_cond_filt$Replicate)
H2AX_quant_allreps_cond_filt$Well <- factor(H2AX_quant_allreps_cond_filt$Well)

# plot H2AX number of spots
# comparisons for statistics
my_comp <- list(c("WT-control", "KO-control"), c("KO-control", "KO-wt-IMPDH2"),
                c("KO-wt-IMPDH2", "KO-3xNLS-IMPDH2"), c("KO-3xNLS-IMPDH2", "KO-NES-IMPDH2"),
                c("KO-control", "KO-3xNLS-IMPDH2"), c("KO-wt-IMPDH2", "KO-NES-IMPDH2"),
                c("KO-control", "KO-NES-IMPDH2"))

# plot number of H2AX spots per nucleus: condition - dots and boxplox
ggplot(H2AX_quant_allreps_cond_filt, aes(x = Condition, y = full.nuclei...Number.of.Spots, fill = Condition)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "", y = "Number of H2AX foci per nucleus", fill = "") +
  scale_fill_manual(values = c("WT-control" = "#3B4E8B", "KO-control" = "#AFC2FF", 
                               "KO-wt-IMPDH2" = "#4C0099",
                               "KO-NES-IMPDH2" = "#A392C0", "KO-3xNLS-IMPDH2" = "#855FFF")) +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")

ggplot(H2AX_quant_allreps_cond_filt, aes(x = Condition, y = full.nuclei...Number.of.Spots, fill = Condition)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "", y = "Number of H2AX foci per nucleus", fill = "") +
  scale_fill_manual(values = c("WT-control" = "#3B4E8B", "KO-control" = "#AFC2FF", 
                               "KO-wt-IMPDH2" = "#4C0099",
                               "KO-NES-IMPDH2" = "#A392C0", "KO-3xNLS-IMPDH2" = "#855FFF")) +
 # stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")

#remove outlier
data_cond_1 <- H2AX_quant_allreps_cond_filt %>%
  group_by(Condition) %>%
  mutate(median_yH2Ax = median(full.nuclei...Number.of.Spots),
         sd_yH2AX = sd(full.nuclei...Number.of.Spots)) %>%
  filter(full.nuclei...Number.of.Spots >= median_yH2Ax - (3 * sd_yH2AX) &
           full.nuclei...Number.of.Spots <= median_yH2Ax + (3 * sd_yH2AX)) %>%
  ungroup()
data_cond_1 <- data_cond_1[!data_cond_1$Condition %in% "KO-NES-IMPDH2",]

my_comp <- list(c("WT-control", "KO-control"), c("KO-control", "KO-wt-IMPDH2"),
                c("KO-wt-IMPDH2", "KO-3xNLS-IMPDH2"),
                c("KO-control", "KO-3xNLS-IMPDH2"))

pdf("/Users/alisc/Downloads/yH2AX_no_outlier.pdf")        
ggplot(data_cond_1, aes(x = Condition, y = full.nuclei...Number.of.Spots, fill = Condition)) +
  geom_boxplot(width = 0.75, fatten = 2) +
  labs(x = "", y = "Number of H2AX foci per nucleus", fill = "") +
  scale_fill_manual(values = c("WT-control" = "#3B4E8B", "KO-control" = "#AFC2FF", 
                               "KO-wt-IMPDH2" = "#4C0099",
                               "KO-NES-IMPDH2" = "#A392C0", "KO-3xNLS-IMPDH2" = "#855FFF")) +
  stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")
dev.off()
