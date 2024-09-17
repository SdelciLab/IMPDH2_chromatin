# libraries
library(tidyverse)
library(ggpubr)


script_name <- basename(rstudioapi::getSourceEditorContext()$path)


setwd("C:/Users/tganez/OneDrive - CRG - Centre de Regulacio Genomica/Toni CRG/P10 - IMPDH2/")



# get data MDAMB231

data_MDA_24h <- read.delim("Fig S3D.txt", skip = 9)

data_MDA_24h <- data_MDA_24h %>%
    mutate(Timepoint = rep("24h", nrow(data_MDA_24h)))

data_label ="imp.parp"


# add conditions
data_MDA_all_cond <- data_MDA_24h %>%
    mutate(Condition = case_when(
        Column == 3 ~ "DMSO",
        Column == 4 ~ "ETO 1uM",
        Column == 5 ~ "ETO 2.5uM",
        Column == 6 ~ "ETO 5uM",
        Column == 7 ~ "ETO 10uM"),
        Well = case_when(
            Row == 5 ~ 1, 
            Row == 6 ~ 2, 
            Row == 7 ~ 3)
    )



# check
table(data_MDA_all_cond$Condition, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Condition, data_MDA_all_cond$Column)

table(data_MDA_all_cond$Well, data_MDA_all_cond$Row)
table(data_MDA_all_cond$Well, data_MDA_all_cond$Column)


# factor
data_MDA_all_cond$Condition <- factor(data_MDA_all_cond$Condition, levels = c("DMSO", "ETO 1uM", "ETO 2.5uM", "ETO 5uM", "ETO 10uM"))
data_MDA_all_cond$Well <- factor(data_MDA_all_cond$Well)





######################################################################
################# insert DAPI analysis for MDA #######################
######################################################################

{
    names(data_MDA_all_cond)[20]='Nuclei_Intensity_Nuc_DAPI'
    names(data_MDA_all_cond)[21]='Nuclei_Ring_Intensity_Nuc_DAPI'
    
}

### Remove background


data_MDA_all_cond$dif_DAPI = data_MDA_all_cond$Nuclei_Intensity_Nuc_DAPI - data_MDA_all_cond$Nuclei_Ring_Intensity_Nuc_DAPI




### Calculate integrated 

data_MDA_all_cond$integrated_DAPI = data_MDA_all_cond$Nuclei.Selected.Selected...Nucleus.Area..µm..*data_MDA_all_cond$dif_DAPI/4250000



### Normalize profiles, adjusted to ETO10h because was the one with recognizable profiles

{
    all_norm=1
    data_MDA_all_cond = data_MDA_all_cond %>%
        mutate(norm_integrated_DAPI=
                   case_when((Condition == "DMSO" & Timepoint =="24h") ~ integrated_DAPI*1.1*all_norm,
                             (Condition == "ETO 1uM" & Timepoint =="24h") ~ integrated_DAPI*1.1*all_norm,
                             (Condition == "ETO 2.5uM" & Timepoint =="24h") ~ integrated_DAPI*1.1*all_norm,
                             (Condition == "ETO 5uM" & Timepoint =="24h") ~ integrated_DAPI*1.1*all_norm,
                             (Condition == "ETO 10uM" & Timepoint =="24h") ~ integrated_DAPI*1.1*all_norm,
                   )
        )
    
    
    ### Set the gates (adjust according to the plots below)
    
    move_all=0
    y_position =1
    gate_seg=0.05
    {
        G1_gate_start =0.55+move_all
        G1_gate_end = 0.95+move_all
        
        G2M_gate_start = 1.25+move_all
        G2M_gate_end = 1.55+move_all
        
        S_gate_start = G1_gate_end+0.01
        S_gate_end = G2M_gate_start-0.01
        
        
        
        custom_gates=c(geom_segment(aes(x = G1_gate_start, xend = G1_gate_end, y = y_position, yend = y_position)),
                       geom_segment(aes(x = G1_gate_start, xend = G1_gate_start, y = y_position+gate_seg, yend = y_position-gate_seg)),
                       geom_segment(aes(x = G1_gate_end, xend = G1_gate_end, y = y_position+gate_seg, yend = y_position-gate_seg)),
                       
                       geom_segment(aes(x = S_gate_start, xend = S_gate_end, y = y_position, yend = y_position)),
                       geom_segment(aes(x = S_gate_start, xend = S_gate_start, y = y_position+gate_seg, yend = y_position-gate_seg)),
                       geom_segment(aes(x = S_gate_end, xend = S_gate_end, y = y_position+gate_seg, yend = y_position-gate_seg)),
                       
                       geom_segment(aes(x = G2M_gate_start, xend = G2M_gate_end, y = y_position, yend = y_position)),
                       geom_segment(aes(x = G2M_gate_start, xend = G2M_gate_start, y = y_position+gate_seg, yend = y_position-gate_seg)),
                       geom_segment(aes(x = G2M_gate_end, xend = G2M_gate_end, y = y_position+gate_seg, yend = y_position-gate_seg)))
        
        custom_gates_vline = geom_vline(xintercept = c(G1_gate_start,G1_gate_end,S_gate_start,S_gate_end,G2M_gate_start,G2M_gate_end),linetype="dotted")
        
    }
    
    
plot_1= ggplot(data_MDA_all_cond, 
               aes(x=norm_integrated_DAPI)) +
    xlim(0.25, 2) +
    geom_density(aes(color =  Timepoint),
                 alpha=0, size=.2, adjust= 0.9
    )+
    
    custom_gates+
    
    ggtitle(paste0("Integrated_DAPI - Eto 10uM\n",script_name))+
    theme_bw()+
                   facet_wrap(Timepoint~Condition,
                              ncol=5)
plot_1
}


############## add the gates data to the main table ###########################

### add the groups to the data table

data_MDA_all_cond = data_MDA_all_cond %>%
    mutate(Phase_Group=                                                 
               case_when((between(norm_integrated_DAPI,G1_gate_start,G1_gate_end)) ~ "G1",
                         (between(norm_integrated_DAPI,S_gate_start,S_gate_end)) ~ "S",
                         (between(norm_integrated_DAPI,G2M_gate_start,G2M_gate_end)) ~ "G2M"))%>%
    filter(Phase_Group %in% c("G1", "S", "G2M")) %>% 
    mutate(Treatment_Phase_Group=paste0(Condition,"\n",Phase_Group))


data_MDA_all_cond$Phase_Group = factor(data_MDA_all_cond$Phase_Group, levels = c("G1", "S", "G2M"))

#####################################################################

# plot facet
my_comp <- list(c("DMSO", "ETO 1uM"), c("ETO 1uM", "ETO 2.5uM"), c("ETO 2.5uM", "ETO 5uM"),
                c("ETO 5uM", "ETO 10uM"))

boxplot = ggboxplot(
    data_MDA_all_cond, x = "Condition", y = "Nuclei.Selected.Selected...Intensity.Nucleus.Alexa.488.Mean", fill = "Condition",
    nrow = 1) +
    stat_compare_means(
        comparisons = my_comp,
        method = "wilcox",
        label = "p.format", size = 2
    )+
    labs(y ="log2(IMPDH2 nuclear mean intensity)") +
    scale_fill_manual(values = c("DMSO" = "#FFE5CC", "ETO 1uM" = "#FFB266",
                                 "ETO 2.5uM" = "#FF9933", "ETO 5uM" = "#FF8000",
                                 "ETO 10uM" = "#CC6600")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust=1))

boxplot +
    facet_wrap(~Phase_Group, ncol = 3)

ggsave(paste("Fig S3D -",data_label,"- log2IMPDH2nuc_MDA_ETO_norm_DMSOfacetstats _ byPhase_transposed_24h.pdf"),  width = 12, height = 6)

