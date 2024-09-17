# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

# get data
filepath = setwd("C:/Users/tganez/OneDrive - CRG - Centre de Regulacio Genomica/Toni CRG/P10 - IMPDH2/")

H2AX_sh_quant <- read.delim("Fig 2F.txt"#, skip = 9
)

# explore
head(H2AX_sh_quant)
colnames(H2AX_sh_quant)

# distribution
table(H2AX_sh_quant$Row, H2AX_sh_quant$Column)

# add condition
H2AX_sh_quant_cond <- H2AX_sh_quant %>%
    mutate(CellLine = case_when(
        Column == 2 | Column == 3 | Column == 4 ~ "MDA-MB-231",
        Column == 6 | Column == 7 | Column == 8 ~ "Cal51"),
        Condition = case_when(
            Column == 2 | Column == 6 ~ "NT",
            Column == 3 | Column == 7 ~ "sh1IMPDH2",
            Column == 4 | Column == 8 ~ "sh2IMPDH2"),
        Replicate = case_when(
            Row == 2 ~ 1,
            Row == 3 ~ 2,
            Row == 4 ~ 3)
    )


# filter interesting data
H2AX_sh_quant_cond_filt <- H2AX_sh_quant_cond %>%
    filter(CellLine %in% c("MDA-MB-231", "Cal51")) %>%
    filter(Condition %in% c("NT", "sh1IMPDH2", "sh2IMPDH2")) %>%
    filter(Replicate %in% c(1,2,3))%>%
    filter(Phase_Group %in% c("G1","S","G2M")) 

# verify
table(H2AX_sh_quant_cond_filt$CellLine, H2AX_sh_quant_cond_filt$Column)
table(H2AX_sh_quant_cond_filt$Condition, H2AX_sh_quant_cond_filt$Column)
table(H2AX_sh_quant_cond_filt$Replicate, H2AX_sh_quant_cond_filt$Row)

# add cell_line + condition
H2AX_sh_quant_cond_filt <- H2AX_sh_quant_cond_filt %>%
    mutate(CellLine_Cond = paste0(CellLine, "_", Condition))

# factors
H2AX_sh_quant_cond_filt$Condition <- factor(H2AX_sh_quant_cond_filt$Condition, levels = 
                                                c("NT", "sh1IMPDH2", "sh2IMPDH2"))
H2AX_sh_quant_cond_filt$CellLine <- factor(H2AX_sh_quant_cond_filt$CellLine, levels = 
                                               c("MDA-MB-231", "Cal51"))
H2AX_sh_quant_cond_filt$CellLine_Cond <- factor(H2AX_sh_quant_cond_filt$CellLine_Cond, levels = 
                                                    c("MDA-MB-231_NT", "MDA-MB-231_sh1IMPDH2", "MDA-MB-231_sh2IMPDH2",
                                                      "Cal51_NT", "Cal51_sh1IMPDH2", "Cal51_sh2IMPDH2"))
H2AX_sh_quant_cond_filt$Replicate <- factor(H2AX_sh_quant_cond_filt$Replicate)

H2AX_sh_quant_cond_filt$Phase_Group <- factor(H2AX_sh_quant_cond_filt$Phase_Group, 
                                               levels = c("G1","S","G2M"))


# comparisons for statistics
my_comp <- list(c("MDA-MB-231_NT", "MDA-MB-231_sh1IMPDH2"), 
                c("MDA-MB-231_NT", "MDA-MB-231_sh2IMPDH2"), 
                c("Cal51_NT", "Cal51_sh1IMPDH2"), 
                c("Cal51_NT", "Cal51_sh2IMPDH2"))

# add condition 
H2AX_sh_quant_cond_filt_2 <- H2AX_sh_quant_cond_filt %>%
    mutate(Condition2 = ifelse(Condition == "NT", "NT", "KD"))

# make the plots

for (line in unique(H2AX_sh_quant_cond_filt$CellLine)) {
    H2AX_sh_quant_cond_filt_3 <- H2AX_sh_quant_cond_filt_2 %>%
        filter(CellLine %in% line)
    
    # comparisons for statistics
    my_comp <- list(c("MDA-MB-231_NT", "MDA-MB-231_sh1IMPDH2"), 
                    c("MDA-MB-231_NT", "MDA-MB-231_sh2IMPDH2"), 
                    c("Cal51_NT", "Cal51_sh1IMPDH2"), 
                    c("Cal51_NT", "Cal51_sh2IMPDH2"))
    
    my_comp=subset(my_comp, grepl(line, my_comp))
    
# final plot: number of spots
plot_1=ggplot(H2AX_sh_quant_cond_filt_3, aes(x = CellLine_Cond, y = full.nuclei...Number.of.Spots, fill = Condition2)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    labs(x = "", y = "Number of H2AX foci per nucleus", fill = "") +
    stat_compare_means(method = "wilcox", size = 2, comparisons = my_comp, tip.length = 0, braket.size = 0.01) +
    scale_fill_manual(name = "", values = c("NT" = "#295C29", "KD" = "#678067"))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position = "none")

plot_1_byPhase = plot_1 +
    facet_wrap(~Phase_Group)
ggsave(paste0("Fig 2F - foci_MDA_", line,"_sh - byPhase.pdf"), device = "pdf", width = 5, height = 5, plot = plot_1_byPhase)


}

# numbers
table(H2AX_sh_quant_cond_filt_2$CellLine_Cond)

