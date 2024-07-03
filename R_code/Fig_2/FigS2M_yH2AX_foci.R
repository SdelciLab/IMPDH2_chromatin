#Figure S2M column of interest: full nuclei - Number of Spots
#WT: B3, B4, B5; KO: C3, C4, C5; KO-WT: D3, D4, D5; KO-CD: E3, E4, E5 file Objects_Population - full nuclei_FigS2M

#load the data
data <- read.delim("/Users/alisc/Downloads/Objects_Population - full nuclei_FigS2M.txt", skip = 9)

data <- data[data$Row != 6,]
data_cond <- data %>%
  mutate(Cell = case_when(
    Row == 2 ~ "WT",
    Row == 3  ~ "KO",
    Row == 4  ~ "KO-WT",
    Row == 5  ~ "KO-CD"
  )
  )
my_comparisons <- list(c("WT", "KO"), c("KO", "KO-WT"), c("KO-WT", "KO-CD"), c("KO", "KO-CD"))
#remove outlier
data_no_outlier <- data_cond %>%
  
  group_by(Cell) %>%
  mutate(median_1 = median(full.nuclei...Number.of.Spots, na.rm= T),
         sd_1 = sd(full.nuclei...Number.of.Spots, na.rm=T)) %>%
  
  filter(full.nuclei...Number.of.Spots >= median_1 - (3 * sd_1) &
           full.nuclei...Number.of.Spots <= median_1 + (3 * sd_1)) %>%
  ungroup() 
ordered <- c("WT", "KO","KO-WT","KO-CD")
data_no_outlier$Cell <- factor(data_no_outlier$Cell, levels = ordered)

pdf("/Users/alisc/Desktop/CRG/final_figures/FigS2M_spots.pdf")
ggplot(data_no_outlier, aes(x = Cell, y = full.nuclei...Number.of.Spots, fill=Cell)) +  
  geom_boxplot(width = 0.75, fatten = 2) + 
  labs(x = "", y = "Number of yH2AX foci per nucleus") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme(legend.position = "none")  + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox")+
  scale_fill_brewer(palette="Accent")
dev.off()

table(data_no_outlier$Cell)
