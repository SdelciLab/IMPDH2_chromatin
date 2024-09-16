# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

# get data
H2AX_RPA_quant <- read.delim("Objects_Population - full nuclei.txt", skip = 9)

# explore
head(H2AX_RPA_quant)
colnames(H2AX_RPA_quant)

# distribution
table(H2AX_RPA_quant$Row, H2AX_RPA_quant$Column)

# add condition
H2AX_RPA_quant_cond <- H2AX_RPA_quant %>%
    mutate(Condition = case_when(
        Column == 2 ~ "WT-control",
        Column == 3 ~ "KO-control",
        Column == 4 ~ "KO-wt-IMPDH2",
        Column == 5 ~ "KO-NES-IMPDH2",
        Column == 6 ~ "KO-3xNLS-IMPDH2"),
        Guanosine = case_when(
            (Row == 2 | Row == 3 | Row == 4) ~ "+Guan",
            (Row == 5 | Row == 6 | Row == 7) ~ "-Guan"),
        Replicate = case_when(
            (Row == 2 | Row == 5) ~ 1,
            (Row == 3 | Row == 6) ~ 2,
            (Row == 4 | Row == 7) ~ 3)
    )

# filter interesting data
H2AX_RPA_quant_cond_filt <- H2AX_RPA_quant_cond %>%
    filter(Condition %in% c("WT-control", "KO-control", "KO-wt-IMPDH2", "KO-NES-IMPDH2", "KO-3xNLS-IMPDH2")) %>%
    filter(Replicate %in% c(1,2,3),
           Guanosine %in% c("-Guan", "+Guan"))

# verify
table(H2AX_RPA_quant_cond_filt$Condition, H2AX_RPA_quant_cond_filt$Column)
table(H2AX_RPA_quant_cond_filt$Replicate, H2AX_RPA_quant_cond_filt$Row)
table(H2AX_RPA_quant_cond_filt$Guanosine, H2AX_RPA_quant_cond_filt$Row)

# factors
H2AX_RPA_quant_cond_filt$Condition <- factor(H2AX_RPA_quant_cond_filt$Condition, levels = 
                                                 c("WT-control", "KO-control", "KO-wt-IMPDH2", "KO-3xNLS-IMPDH2",
                                                   "KO-NES-IMPDH2"))
H2AX_RPA_quant_cond_filt$Replicate <- factor(H2AX_RPA_quant_cond_filt$Replicate)
H2AX_RPA_quant_cond_filt$Guanosine <- factor(H2AX_RPA_quant_cond_filt$Guanosine, levels = 
                                                 c("-Guan", "+Guan"))

# select data WT_KO
H2AX_RPA_quant_cond_WT_KO <- H2AX_RPA_quant_cond_filt %>%
    filter((Condition == "WT-control" & Guanosine == "-Guan") | Condition == "KO-control") %>%
    mutate(Cond_Guan = paste0(Condition, "_", Guanosine))
table(H2AX_RPA_quant_cond_WT_KO$Condition, H2AX_RPA_quant_cond_WT_KO$Guanosine)
table(H2AX_RPA_quant_cond_WT_KO$Cond_Guan)

# factor
H2AX_RPA_quant_cond_WT_KO$Cond_Guan <- factor(H2AX_RPA_quant_cond_WT_KO$Cond_Guan, levels = 
                                                  c("WT-control_-Guan", "KO-control_-Guan", "KO-control_+Guan"))

# get comparisons
my_comp <- list(c("WT-control_-Guan", "KO-control_-Guan"), c("KO-control_-Guan", "KO-control_+Guan"),
                c("WT-control_-Guan", "KO-control_+Guan"))

# plot area nucleus: full.nuclei...areainum2n.Area..µm..
ggplot(H2AX_RPA_quant_cond_WT_KO, aes(x = Cond_Guan, y = log2(full.nuclei...areainum2n.Area..µm..), 
                                      fill = Cond_Guan)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "", y = "log2(Area nuclei (um2))", fill = "") +
    scale_fill_manual(name = "", labels = c("WT", "KO"),
                      values = c("WT-control_-Guan" = "#3B4E8B", "KO-control_-Guan" = "#AFC2FF",
                                 "KO-control_+Guan" = "#AFC2FF"))+
    stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")
ggsave("plots/log2areanuclei_WT_KO.pdf", device = "pdf", width = 3.5, height = 4)

# numbers
table(H2AX_RPA_quant_cond_WT_KO$Cond_Guan)