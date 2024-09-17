


library(dplyr)
library(sp)
library(ggplot2)
library(ggpubr)
library(magick)
library(ggforce)

script_name <- basename(rstudioapi::getSourceEditorContext()$path)

########################################################################
###### prepare data ####################################################
########################################################################


filepath = setwd("C:/Users/tganez/OneDrive - CRG - Centre de Regulacio Genomica/Toni CRG/P10 - IMPDH2/")


# load the files and remove columns that are all NA
filename_object = 
  read.csv("Fig 2D.txt", header = TRUE, sep="\t",skip=9) %>%
  select(where(function(x) any(!is.na(x)))) 


nuclei_data = filename_object %>%
    mutate(Antibody=
             case_when(
                       (between(Column,7,9) & Row==7) ~ "gH2AX"
             ),
           Cell_line=
             case_when((between(Column,1,12) & (between(Row,1,8))) ~ "U2OS_FUCCI"
             ),
           Biological_rep=
             case_when((Column %in% c(1,4,7,10)) ~ "Rep1",
                       (Column %in% c(2,5,8,11)) ~ "Rep2",
                       (Column %in% c(3,6,9,12)) ~ "Rep3"
             ),
    )

 
# get the interesting columns
nuclei_data= nuclei_data[,c(16:38)]

{
names(nuclei_data)[1]='Object n sel sel'
names(nuclei_data)[2]='Nuc_Area'
names(nuclei_data)[3]='Nuc_Roundness'
names(nuclei_data)[4]='Cyt_Area'
names(nuclei_data)[5]='Cyt_Roundness'
names(nuclei_data)[6]='Cell_Area'
names(nuclei_data)[7]='Cell_Roundness'
names(nuclei_data)[8]='NAdivCAxCR'
names(nuclei_data)[9]='nuc_int_488'
names(nuclei_data)[10]='nuc_int_Turq'
names(nuclei_data)[11]='nuc_int_546'
names(nuclei_data)[12]='nuc_int_647'
names(nuclei_data)[13]='cyt_int_488'
names(nuclei_data)[14]='cyt_int_Turq'
names(nuclei_data)[15]='cyt_int_546'
names(nuclei_data)[16]='cyt_int_647'
names(nuclei_data)[17]='nuc_intR_488'
names(nuclei_data)[18]='nuc_intR_Turq'
names(nuclei_data)[19]='nuc_intR_546'
names(nuclei_data)[20]='nuc_intR_647'
}
            


########################################
##### adjust FUCCI plot ################    
########################################

{
  ### substract green signal from around the nucleus 
  nuclei_data$dif488=nuclei_data$nuc_int_488 - nuclei_data$nuc_intR_488
  nuclei_data$difTurq=nuclei_data$nuc_int_Turq - nuclei_data$nuc_intR_Turq
  nuclei_data$dif546=nuclei_data$nuc_int_546 - nuclei_data$nuc_intR_546
  nuclei_data$nuc_int_647_intgr = nuclei_data$nuc_int_647*nuclei_data$Nuc_Area
  
  ### compensate bleed-through
  
  m_var = 0.89
  c_var = -56.59
  nuclei_data$dif488_int_btcomp = nuclei_data$dif488-(m_var*nuclei_data$difTurq+c_var) # removing the constant the fucci profile/yellow makes more sense
  

  
  move_Turq = -0.25
  stretch_Turq = 1.15
  move_488= -1.25
  stretch_488 = 1.65

  ### center the data in the gates
  nuclei_data$difTurq_mod = nuclei_data$difTurq # "+-" for moving whole set, "*" for stretching and squishing
  nuclei_data$dif488_int_btcomp_mod = nuclei_data$dif488_int_btcomp
  
#####
  
  nuclei_data$log10_difTurq_mod = (log10(nuclei_data$difTurq_mod) + move_Turq) * stretch_Turq
  nuclei_data$log10_dif488_int_btcomp_mod = (log10(nuclei_data$dif488_int_btcomp_mod) + move_488) * stretch_488
  
  
  ### make the gates
  
  move_gates_x = .75
  stretch_gates_x = 0.9
  move_gates_y= -0.55
  stretch_gates_y = 1.1
  
  
  My = (c(2.675, 2.675, 3.800, 3.800) + move_gates_y) * stretch_gates_y
  Mx = (c(0.000, 1.619, 1.619, 0.000) + move_gates_x) * stretch_gates_x
  
  G2My = (c(2.675, 2.676, 3.376, 3.900, 3.900, 3.800) + move_gates_y) * stretch_gates_y
  G2Mx = (c(1.620, 1.625, 3.000, 3.000, 2.375, 1.620) + move_gates_x) * stretch_gates_x
  
  Sy =   (c(2.001, 2.001, 2.226, 2.506, 3.375, 2.675) + move_gates_y) * stretch_gates_y
  Sx =   (c(1.625, 1.875, 2.615, 3.000, 3.000, 1.625) + move_gates_x) * stretch_gates_x
  
  G1y =  (c(1.250, 1.125, 1.175, 2.225, 2.000, 2.000) + move_gates_y) * stretch_gates_y
  G1x =  (c(0.500, 1.875, 2.815, 2.615, 1.875, 0.500) + move_gates_x) * stretch_gates_x
  
  G0y =  (c(1.175, 1.675, 2.226, 2.875, 2.225) + move_gates_y) * stretch_gates_y
  G0x =  (c(2.816, 3.625, 3.750, 3.625, 2.616) + move_gates_x) * stretch_gates_x
  
  
  ###################################
  ###### make the test plot #########
  ###################################
  
  
  nuclei_data_filt=nuclei_data %>%
    filter(Antibody %in% c("gH2AX"))%>%
    filter(Cell_line %in% c("U2OS_FUCCI"))
  
  voro_nuc546plot_GT = ggplot(nuclei_data_filt %>%
                                mutate(avg_ant =mean(nuc_int_546) ) %>%
                                mutate(nuc_int_546_norm = nuc_int_546/avg_ant) , aes(x=log10_difTurq_mod, y=log10_dif488_int_btcomp_mod)) +
    geom_voronoi_tile(aes(fill = log2(nuc_int_546_norm)), max.radius = 0.05) +
    geom_polygon(data = data.frame(G1x, G1y),aes(x = G1x, y = G1y), colour = "black", fill=NA, linetype = "longdash")+
    geom_polygon(data = data.frame(Sx, Sy),aes(x = Sx, y = Sy), colour = "black", fill=NA ,linetype = "longdash")+
    geom_polygon(data = data.frame(G2Mx, G2My),aes(x = G2Mx, y = G2My), colour = "black", fill=NA ,linetype = "longdash")+
    geom_polygon(data = data.frame(G0x, G0y),aes(x = G0x, y = G0y), colour = "black", fill=NA ,linetype = "longdash")+
    geom_polygon(data = data.frame(Mx, My),aes(x = Mx, y = My), colour = "black", fill=NA ,linetype = "longdash")+
    theme_bw()+
    scale_fill_gradient2(low = "blue",
                         mid = "khaki",
                         high = "red",
                         na.value = NA,#limits=c(-2,2))+
                         limits= c(-1,1), oob = scales::squish)+
    xlab("Turquoise2")+
    ylab("Clover")+
    coord_cartesian(ylim = c(0,4),
                    xlim = c(0,4))+
    facet_wrap(~Antibody)+
    ggtitle(paste("- Voronoi 546"))
  
  print(voro_nuc546plot_GT)
}

# save gates coordinates in a dataframe for points in polygon
  allx = c(G0x,G1x,Sx,G2Mx,Mx)
  ally = c(G0y,G1y,Sy,G2My,My)
  gatenames = c(rep("G0",length(G0x)),rep("G1",length(G1x)),rep("S",length(Sx)),rep("G2M",length(G2Mx)),rep("M",length(Mx)))
  gates = data.frame(allx,ally,gatenames)

### Get the points that fall in each gate
  G0gate = point.in.polygon(nuclei_data$log10_difTurq_mod,nuclei_data$log10_dif488_int_btcomp_mod,G0x,G0y)
  G1gate = point.in.polygon(nuclei_data$log10_difTurq_mod,nuclei_data$log10_dif488_int_btcomp_mod,G1x,G1y)
  Sgate = point.in.polygon(nuclei_data$log10_difTurq_mod,nuclei_data$log10_dif488_int_btcomp_mod,Sx,Sy)
  G2Mgate = point.in.polygon(nuclei_data$log10_difTurq_mod,nuclei_data$log10_dif488_int_btcomp_mod,G2Mx,G2My)
  Mgate = point.in.polygon(nuclei_data$log10_difTurq_mod,nuclei_data$log10_dif488_int_btcomp_mod,Mx,My)


### Add columns to say in which gate they are present (1 present, 0 not present)
  nuclei_data$G0gate = G0gate
  nuclei_data$G1gate = G1gate
  nuclei_data$Sgate = Sgate
  nuclei_data$G2Mgate = G2Mgate
  nuclei_data$Mgate = Mgate
  
  nuclei_data = nuclei_data %>%
    mutate(Phase_Group=
             case_when((G0gate==1)  ~ "G0",
                       (G1gate==1)  ~ "G1",
                       (Sgate==1)  ~ "S",
                       (G2Mgate==1)  ~ "G2M",
                       (Mgate==1)  ~ "M"
                       )
    )


  
  ###################################
  ###### make the final plots #######
  ###################################
 
targets =c("gH2AX")
 
  nuclei_data_filt=nuclei_data %>%
    filter(Antibody %in% c(targets))%>%
    filter(Cell_line %in% c("U2OS_FUCCI"))
  
    # force the order of the factors to be the same as the order in the vector
    nuclei_data_filt$Antibody = factor(nuclei_data_filt$Antibody,levels=c(targets))
    
 
### 647 Integrated   
    
    targets =c("gH2AX")
    
    plot_list=c()
    
    for (target in targets) {
        nuclei_data_filt=nuclei_data %>%
        filter(Antibody %in% c(target))%>%
        filter(Cell_line %in% c("U2OS_FUCCI"))
      
      # force the order of the factors to be the same as the order in the vector
      nuclei_data_filt$Antibody = factor(nuclei_data_filt$Antibody,levels=c(target,"neg"))
      
      voro_nuc647plot_GT = ggplot(nuclei_data_filt %>%
                                    mutate(avg_ant =mean(nuc_int_647_intgr) ) %>%
                                    mutate(nuc_int_647_norm = nuc_int_647_intgr/avg_ant) , aes(x=log10_difTurq_mod, y=log10_dif488_int_btcomp_mod)) +
        geom_voronoi_tile(aes(fill = log2(nuc_int_647_norm)), max.radius = 0.05) +
        
        geom_polygon(data = data.frame(G1x, G1y),aes(x = G1x, y = G1y), colour = "black", fill=NA, linetype = "longdash")+
        geom_polygon(data = data.frame(Sx, Sy),aes(x = Sx, y = Sy), colour = "black", fill=NA ,linetype = "longdash")+
        geom_polygon(data = data.frame(G2Mx, G2My),aes(x = G2Mx, y = G2My), colour = "black", fill=NA ,linetype = "longdash")+
        geom_polygon(data = data.frame(G0x, G0y),aes(x = G0x, y = G0y), colour = "black", fill=NA ,linetype = "longdash")+
        geom_polygon(data = data.frame(Mx, My),aes(x = Mx, y = My), colour = "black", fill=NA ,linetype = "longdash")+
        theme_bw()+
        scale_fill_gradient2(low = "blue",
                             mid = "khaki",
                             high = "red",
                             na.value = NA,
                             limits= c(-1,1), 
                             oob = scales::squish,
                             name = "norm\n647")+
        xlab("Turquoise2")+
        ylab("Clover")+
        coord_cartesian(ylim = c(0,4.6),
                        xlim = c(0,4.6))+
        facet_wrap(~Antibody)+
        ggtitle(paste0(script_name,"\n",target," - Voronoi 647 integrated"))
      
      print(voro_nuc647plot_GT)
    }
    

    ##############################################################
    ################ boxplot #####################################
    ##############################################################

    filtered_data = nuclei_data %>%
      filter(Antibody %in% c(target))%>%
      filter(!is.na(Phase_Group)) %>% 
      group_by(Antibody) %>% 
      mutate(avg_ant =mean(nuc_int_647_intgr)) %>%
      mutate(nuc_int_647_norm = nuc_int_647_intgr/avg_ant) %>% 
      filter(Phase_Group %in% c("G1", "S", "G2M")) %>%
      na.omit()
    
    filtered_data$Phase_Group= factor(filtered_data$Phase_Group,levels = c("G1","S","G2M"))
    
    # Count the number of items for each Phase Group
    table(filtered_data$Phase_Group)
    
    my_comp <- list(c("G1", "S"), 
                    c("S", "G2M"), 
                    c("G1", "G2M")
    )
    
    set.seed(123)
    boxplot <- ggplot(na.omit(filtered_data %>% 
                                    filter(Phase_Group %in% c("G1", "S", "G2M"))
    ), aes(x = Phase_Group, y = log2(nuc_int_647_norm))) +
      geom_jitter(color="grey", shape = 16, alpha = 0.50) +
      geom_boxplot(fill = NA, lwd = 0.75) +
      
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black"),
            strip.background = element_rect(colour = "black"),
            plot.title = element_text(size = 10)) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = my_comp,
        label = "p.format", size = 2
      )+
      ggtitle(paste0(script_name, "\nnuc_int_647_intgr - gH2AX"))
    
    print(boxplot)
    ggsave(paste("Fig 2D - nuclear gH2AX integrated intensity.pdf"),
           width = 12 + 0, height = 12 + 0, dpi = 150, units = "cm")
    
