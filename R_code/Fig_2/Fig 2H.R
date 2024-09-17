# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

filepath = setwd("C:/Users/tganez/OneDrive - CRG - Centre de Regulacio Genomica/Toni CRG/P10 - IMPDH2/")
# get data
H2AX_MPA_quant <- read.delim("Fig 2H.txt")

# explore
head(H2AX_MPA_quant)
colnames(H2AX_MPA_quant)

# distribution
table(H2AX_MPA_quant$Row, H2AX_MPA_quant$Column)

# add condition
H2AX_MPA_quant_cond <- H2AX_MPA_quant %>%
    mutate(Condition = case_when(
        Column == 2 ~ "MPA 0 uM",
        Column == 3 ~ "MPA 2.5 uM",
        Column == 4 ~ "MPA 7.5 uM",
        Column == 5 ~ "MPA 15 uM",
        Column == 7 ~ "DMSO",
        Column == 8 ~ "ETO"),
        Replicate = case_when(
            Row == 2 ~ 1,
            Row == 3 ~ 2,
            Row == 4 ~ 3)
    )

# filter interesting data
H2AX_MPA_quant_cond_filt <- H2AX_MPA_quant_cond %>%
    filter(Condition %in% c("MPA 0 uM", "MPA 2.5 uM", "MPA 7.5 uM", "MPA 15 uM", "DMSO", "ETO")) %>%
    filter(Replicate %in% c(1,2,3))  %>%
    filter(Phase_Group %in% c("G1","S","G2M")) 

# verify
table(H2AX_MPA_quant_cond_filt$Condition, H2AX_MPA_quant_cond_filt$Column)
table(H2AX_MPA_quant_cond_filt$Replicate, H2AX_MPA_quant_cond_filt$Row)

# factors
H2AX_MPA_quant_cond_filt$Condition <- factor(H2AX_MPA_quant_cond_filt$Condition, 
                                             levels = c("MPA 0 uM", "MPA 2.5 uM", "MPA 7.5 uM", "MPA 15 uM", "DMSO", "ETO"))
H2AX_MPA_quant_cond_filt$Replicate <- factor(H2AX_MPA_quant_cond_filt$Replicate)
H2AX_MPA_quant_cond_filt$Phase_Group <- factor(H2AX_MPA_quant_cond_filt$Phase_Group, 
                                             levels = c("G1","S","G2M"))


# remove 3rd rep DMSO because it was unusually high
H2AX_MPA_quant_cond_filt_noout <- H2AX_MPA_quant_cond_filt %>%
    filter(!(Condition == "DMSO" & Replicate == 3))

# comparisons for statistics
my_comp <- list(c("MPA 0 uM", "MPA 2.5 uM"), 
                c("MPA 2.5 uM", "MPA 7.5 uM"), 
                c("MPA 7.5 uM", "MPA 15 uM"),
                c("MPA 0 uM", "MPA 7.5 uM"), 
                c("MPA 2.5 uM", "MPA 15 uM"), 
                c("MPA 0 uM", "MPA 15 uM"),
                c("DMSO", "ETO"))

# plot number of H2AX foci
ggplot_1=ggplot(H2AX_MPA_quant_cond_filt_noout, aes(x = Condition, y = full.nuclei...Number.of.Spots, 
                                           fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "", y = "Number of H2AX foci per nucleus", fill = "") +
    scale_fill_manual(values = c("MPA 0 uM" = "#FFE5CC", "MPA 2.5 uM" = "#FFB266", 
                                 "MPA 7.5 uM" = "#FF9933", "MPA 15 uM" = "#FF8000", 
                                 "DMSO" = "#FFE5CC", "ETO" = "#CC6600")) +
    stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")

ggplot_1_byPhase=ggplot_1+
    facet_wrap(~Phase_Group)
ggsave("Fig 2H - foci_MPA-byPhase.pdf", device = "pdf", width = 5, height = 5, ggplot_1_byPhase) #ggsave("plots/nfoci_MPA.pdf", device = "pdf", width = 5, height = 5)

# number of cells
table(H2AX_MPA_quant_cond_filt_noout$Condition)
