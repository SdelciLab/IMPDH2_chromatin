


library(dplyr)
library(sp)
library(ggplot2)
library(ggpubr)
library(gganimate)
library(magick)
library(ggforce)
library(ggh4x)

# Get the path of the currently active script
script_path <- rstudioapi::getSourceEditorContext()$path

# Extract the script name from the path
script_name <- basename(script_path)

### Prepare data

filepath = setwd("C:/Users/tganez/OneDrive - CRG - Centre de Regulacio Genomica/Toni CRG/P10 - IMPDH2/")

{
  # load the files and remove columns that are all NA
  filename_object = 
    read.csv("Fig 2A 2B S2A.txt", header = TRUE, sep="\t",skip=9)
  

  fuccidata = filename_object
  fuccidata = fuccidata[,c(1,2,4:7,20:39)]
  
  # change the name of the columns
  {
    names(fuccidata)[1]='Row'
    names(fuccidata)[2]='Column'
    names(fuccidata)[3]='Field'
    names(fuccidata)[4]='Object n'
    names(fuccidata)[5]='posX'
    names(fuccidata)[6]='posY'
    names(fuccidata)[7]='Object n sel sel'
    names(fuccidata)[8]='Nuc_Area'
    names(fuccidata)[9]='Nuc_Roundness'
    names(fuccidata)[10]='Cyt_Area'
    names(fuccidata)[11]='Cyt_Roundness'
    names(fuccidata)[12]='Cell_Area'
    names(fuccidata)[13]='Cell_Roundness'
    names(fuccidata)[14]='NAdivCAxCR'
    names(fuccidata)[15]='nuc_int_488'
    names(fuccidata)[16]='nuc_int_Turq'
    names(fuccidata)[17]='nuc_int_546'
    names(fuccidata)[18]='nuc_int_647'
    names(fuccidata)[19]='cyt_int_488'
    names(fuccidata)[20]='cyt_int_Turq'
    names(fuccidata)[21]='cyt_int_546'
    names(fuccidata)[22]='cyt_int_647'
    names(fuccidata)[23]='nuc_intR_488'
    names(fuccidata)[24]='nuc_intR_Turq'
    names(fuccidata)[25]='nuc_intR_546'
    names(fuccidata)[26]='nuc_intR_647'
  }
} 

{
### substract green signal from around the nucleus 
    fuccidata$dif488 = fuccidata$nuc_int_488 - fuccidata$nuc_intR_488
    fuccidata$difTurq = fuccidata$nuc_int_Turq - fuccidata$nuc_intR_Turq
    fuccidata$dif546 = fuccidata$nuc_int_546 - fuccidata$nuc_intR_546

### compensate bleed-through
  fuccidata$dif488_int_btcomp = fuccidata$dif488-(0.56*fuccidata$difTurq-56.59) # removing the constant the fucci profile/yellow makes more sense
  
### calculate integrated intensities
  fuccidata$nuc_int_647_intgr = fuccidata$nuc_int_647*fuccidata$Nuc_Area
  fuccidata$cyt_int_647_intgr = fuccidata$cyt_int_647*fuccidata$Cyt_Area
  
### center the data in the gates
  fuccidata$difTurq_mod = fuccidata$difTurq # "+-" for moving whole set, "*" for stretching and squishing
  fuccidata$dif488_int_btcomp_mod = fuccidata$dif488_int_btcomp
  
  
  move_Turq = -1.15 # more negative moves data to the left
  stretch_Turq = 1.25
  move_488= -.75 # higher moves up
  stretch_488 = 1.2 #higher spreads data
}
  
  fuccidata$log10_difTurq_mod = (log10(fuccidata$difTurq_mod) + move_Turq) * stretch_Turq
  fuccidata$log10_dif488_int_btcomp_mod = (log10(fuccidata$dif488_int_btcomp_mod) + move_488) * stretch_488
  
  
  ### make the gates
  move_gates_x = .55
  stretch_gates_x = .85
  move_gates_y= -0.55
  stretch_gates_y = 1.1
  
  
  My = (c(2.675, 2.675, 3.800, 3.800) + move_gates_y) * stretch_gates_y
  Mx = (c(0.000, 1.619, 1.619, 0.000) + move_gates_x) * stretch_gates_x
  
  G2My = (c(2.675, 2.676, 3.376, 3.900, 3.900, 3.800) + move_gates_y) * stretch_gates_y
  G2Mx = (c(1.620, 1.625, 3.000, 3.000, 2.375, 1.620) + move_gates_x) * stretch_gates_x
  
  Sy =  (c(2.001, 2.001, 2.226, 2.506, 3.375, 2.675) + move_gates_y) * stretch_gates_y
  Sx =  (c(1.625, 1.875, 2.615, 3.000, 3.000, 1.625) + move_gates_x) * stretch_gates_x
 
  G1y = (c(1.250, 1.125, 1.275, 2.225, 2.000, 2.000) + move_gates_y) * stretch_gates_y
  G1x = (c(0.500, 1.875, 2.615, 2.615, 1.875, 0.500) + move_gates_x) * stretch_gates_x
  
  G0y = (c(1.275, 1.675, 2.226, 2.875, 2.225) + move_gates_y) * stretch_gates_y
  G0x = (c(2.616, 3.625, 3.750, 3.625, 2.616) + move_gates_x) * stretch_gates_x
 
  
  # Save gates coordinates in a dataframe for points in polygon
  allx = c(G0x,G1x,Sx,G2Mx,Mx)
  ally = c(G0y,G1y,Sy,G2My,My)
  gatenames = c(rep("G0",length(G0x)),rep("G1",length(G1x)),rep("S",length(Sx)),rep("G2M",length(G2Mx)),rep("M",length(Mx)))
  gates = data.frame(allx,ally,gatenames)
  
  # Get the points that fall in each gate
  G0gate = point.in.polygon(fuccidata$log10_difTurq_mod,fuccidata$log10_dif488_int_btcomp_mod,G0x,G0y)
  G1gate = point.in.polygon(fuccidata$log10_difTurq_mod,fuccidata$log10_dif488_int_btcomp_mod,G1x,G1y)
  Sgate = point.in.polygon(fuccidata$log10_difTurq_mod,fuccidata$log10_dif488_int_btcomp_mod,Sx,Sy)
  G2Mgate = point.in.polygon(fuccidata$log10_difTurq_mod,fuccidata$log10_dif488_int_btcomp_mod,G2Mx,G2My)
  Mgate = point.in.polygon(fuccidata$log10_difTurq_mod,fuccidata$log10_dif488_int_btcomp_mod,Mx,My)
  
  
  # Gates where they are present (1 present, 0 not present)
  fuccidata$G0gate = G0gate
  fuccidata$G1gate = G1gate
  fuccidata$Sgate = Sgate
  fuccidata$G2Mgate = G2Mgate
  fuccidata$Mgate = Mgate
  
  fuccidata = fuccidata %>%
    mutate(Phase_Group=
             case_when((G0gate==1)  ~ "G0",
                       (G1gate==1)  ~ "G1",
                       (Sgate==1)  ~ "S",
                       (G2Mgate==1)  ~ "G2M",
                       (Mgate==1)  ~ "M"
             )
    )

  
topbot_percent = 10

{
fuccidata_filt_top_c=fuccidata %>%
  arrange(cyt_int_647_intgr)%>%
  slice_tail(prop = topbot_percent/100)


fuccidata_filt_bottom_c=fuccidata %>%
  arrange(cyt_int_647_intgr)%>%
  slice_head(prop = topbot_percent/100)

fuccidata_filt_top_c$topbot ="top"
fuccidata_filt_bottom_c$topbot ="bot"
fuccidata_filt_topbotc = rbind(fuccidata_filt_top_c,fuccidata_filt_bottom_c)

numbers=c(nrow(fuccidata),nrow(fuccidata_filt_top_c),nrow(fuccidata_filt_bottom_c))
n=data.frame(numbers)
n$topbotn=c(paste0("total "),  
           paste0("top ",topbot_percent,"% "),
           paste0("bottom ",topbot_percent,"% "))
n$labels=paste(n$topbotn,"n =",n$number)
n$pos_x = c(0.7, 0.7, 0.7)
n$pos_y = c(0.5, 0.3, 0.1)
}


### Density plots
{
  cyt647plot_GT_topbot_density=ggplot(fuccidata_filt_topbotc, aes(x=log10_difTurq_mod, y=log10_dif488_int_btcomp_mod)) +
    stat_density_2d(geom = "polygon",
                    aes(colour = topbot, fill=topbot), 
                    bins = 6,
                    contour_var = 'ndensity',
                    contour = TRUE,
                    alpha = 0.1) + 
    scale_fill_manual(values=c("top"="red", "bot"="blue")) + 
    scale_color_manual(values=c("top"="red", "bot"="blue")) +
    geom_polygon(data = data.frame(G1x, G1y),aes(x = G1x, y = G1y), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(Sx, Sy),aes(x = Sx, y = Sy), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(G2Mx, G2My),aes(x = G2Mx, y = G2My), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(G0x, G0y),aes(x = G0x, y = G0y), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(Mx, My),aes(x = Mx, y = My), colour = "black", fill=NA,linetype = "longdash")+
    theme_classic()+
    xlab("Turquoise2")+
    ylab("Clover")+
    coord_cartesian(ylim = c(0,4),
                    xlim = c(0,4))+
    geom_text(data=n,
              mapping= aes(x=pos_x,
                           y=pos_y,
                           hjust = "outward",
                           label=labels))+
    facet_wrap(~topbot)+
    ggtitle(paste(script_name,"\nU2OS - Topbot",topbot_percent,"% cyt647plot_GT_topbot_density"))
  
 print(cyt647plot_GT_topbot_density) 
}


{
fuccidata_filt_top_n=fuccidata %>%
  arrange(nuc_int_647_intgr)%>%
  slice_tail(prop = topbot_percent/100)


fuccidata_filt_bottom_n=fuccidata %>%
  arrange(nuc_int_647_intgr)%>%
  slice_head(prop = topbot_percent/100)

fuccidata_filt_top_n$topbot ="top"
fuccidata_filt_bottom_n$topbot ="bot"
fuccidata_filt_topbotn = rbind(fuccidata_filt_top_n,fuccidata_filt_bottom_n)


numbers=c(nrow(fuccidata),nrow(fuccidata_filt_top_n),nrow(fuccidata_filt_bottom_n))
n=data.frame(numbers)
n$topbotn=c(paste0("total "),
           paste0("top ",topbot_percent,"% "),
           paste0("bottom ",topbot_percent,"% "))
n$labels=paste(n$topbotn,"n =",n$number)
n$pos_x = c(0.7, 0.7, 0.7)
n$pos_y = c(0.5, 0.3, 0.1)
  
{
  nuc647plot_GT_topbot_density=ggplot(fuccidata_filt_topbotn
                                       , aes(x=log10_difTurq_mod, y=log10_dif488_int_btcomp_mod)) +
    stat_density_2d(geom = "polygon",
                    aes(colour = topbot, fill=topbot),
                    bins = 6,
                    contour_var = 'ndensity',
                    contour = TRUE,
                    alpha = 0.1) +
    scale_fill_manual(values=c("top"="red", "bot"="blue")) +
    scale_color_manual(values=c("top"="red", "bot"="blue")) +
    geom_polygon(data = data.frame(G1x, G1y),aes(x = G1x, y = G1y), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(Sx, Sy),aes(x = Sx, y = Sy), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(G2Mx, G2My),aes(x = G2Mx, y = G2My), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(G0x, G0y),aes(x = G0x, y = G0y), colour = "black", fill=NA,linetype = "longdash")+
    geom_polygon(data = data.frame(Mx, My),aes(x = Mx, y = My), colour = "black", fill=NA,linetype = "longdash")+
    theme_classic()+
    xlab("Turquoise2")+
    ylab("Clover")+
    coord_cartesian(ylim = c(0,4),
                    xlim = c(0,4))+
    facet_wrap(~topbot)+
    geom_text(data=n,
              mapping= aes(x=pos_x,
                           y=pos_y,
                           hjust = "outward",
                           label=labels))+
    ggtitle(paste(script_name,"\nU2OS - Topbot",topbot_percent,"% nuc647plot_GT_topbot_density"))
  
 print(nuc647plot_GT_topbot_density)
}
 pdfname = paste("Fig S2A - U2OS - Topbot",topbot_percent," cyt647plot_GT_topbot_density.pdf")
  ggsave(pdfname, plot = cyt647plot_GT_topbot_density, width = 40, height = 15, dpi = 150, units = "cm")

 pdfname = paste("Fig 2A - U2OS - Topbot",topbot_percent," nuc647plot_GT_topbot_density.pdf")
  ggsave(pdfname, plot = nuc647plot_GT_topbot_density, width = 40, height = 15, dpi = 150, units = "cm")
}




### Boxplots

filtered_data <- fuccidata %>%
  filter(Phase_Group %in% c("G1", "S", "G2M")) %>%
  na.omit()

filtered_data$Phase_Group = factor(filtered_data$Phase_Group,levels =c("G1","S","G2M"))

# Count the number of items for each Phase Group
table(filtered_data$Phase_Group)

my_comp <- list(c("G1", "S"), 
                c("S", "G2M"), 
                c("G1", "G2M")
                )

    set.seed(123)
    boxplot_all <- ggplot(na.omit(filtered_data), aes(x = Phase_Group, y = log2(nuc_int_647_intgr))) +
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
      ggtitle(paste0(script_name, "\nnuc_int_647_intgr - All"))
    
    print(boxplot_all)
      ggsave(paste("Fig 2B - nuclear IMPDH2 integrated intensity.pdf"),
             width = 12 + 0, height = 12 + 0, dpi = 150, units = "cm")




