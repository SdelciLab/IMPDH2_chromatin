library(DESeq2)
library(ggplot2)
library(tidyr)
library(tidyverse)

#load the DeSeq data
load("/Users/alisc/Desktop/CRG/DeSeq R/DDS_file.RData")
count_data <- counts(dds)
rld <- rlog(dds, blind = TRUE) #normalize the counts first
rld_counts <- assay(rld)

SIRT6_counts <- rld_counts[rownames(count_data) == "ENSG00000077463", ]
SIRT6_counts2 <- as.data.frame(SIRT6_counts)
SIRT6_counts2$condition <- rownames(SIRT6_counts2)


SIRT6_counts2 <- SIRT6_counts2 %>%
  mutate(condition = sub("[0-9]+$", "", condition))
SIRT6_counts2 <- SIRT6_counts2[!SIRT6_counts2$condition %in% c("K0NES"),]
SIRT6_counts2_mean <- SIRT6_counts2 %>%
  group_by(condition) %>%
  summarise(mean(SIRT6_counts), sd(SIRT6_counts))
colnames(SIRT6_counts2_mean) <- c("condition", "Mean", "Std")

ggplot(SIRT6_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean SIRT6 Counts with Standard Deviation",
       y = "Mean normalised SIRT6 Counts")

order2 <- c("WT", "K0WT", "K0NLS", "K0SIG", "K0N0G")
SIRT6_counts2_mean$condition <- factor(SIRT6_counts2_mean$condition, levels = order2)
ggplot(SIRT6_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean SIRT6 Counts with Standard Deviation",
       y = "Mean normalised SIRT6 Counts")

#t test
wt_data <- SIRT6_counts2 %>% filter(condition == "WT") %>% select(SIRT6_counts)
kowt_data <- SIRT6_counts2 %>% filter(condition == "K0WT") %>% select(SIRT6_counts)
konls_data <- SIRT6_counts2 %>% filter(condition == "K0NLS") %>% select(SIRT6_counts)
konog_data <- SIRT6_counts2 %>% filter(condition == "K0N0G") %>% select(SIRT6_counts)
kosig_data <- SIRT6_counts2 %>% filter(condition == "K0SIG") %>% select(SIRT6_counts)

t_test_WT_KOWT <- t.test(wt_data$SIRT6_counts, kowt_data$SIRT6_counts, paired = FALSE)
t_test_KOWT_KONLS <- t.test(konls_data$SIRT6_counts, kowt_data$SIRT6_counts, paired = FALSE)
t_test_KOWT_KOSIG <- t.test(kosig_data$SIRT6_counts, kowt_data$SIRT6_counts, paired = FALSE)
t_test_KOWT_KONOG <- t.test(konog_data$SIRT6_counts, kowt_data$SIRT6_counts, paired = FALSE)


ggplot(SIRT6_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean SIRT6 Counts with Standard Deviation",
       y = "Mean normalised SIRT6 Counts") + geom_segment(aes(x = 1, xend = 2, y = 10.01, yend =10.01)) + 
  geom_segment(aes(x = 2, xend = 3, y = 10.015, yend = 10.015)) + geom_segment(aes(x = 2, xend = 4, y = 10.1, yend =10.1)) + 
  geom_segment(aes(x = 2, xend = 5, y = 10.2, yend = 10.2)) + 
  annotate("text", x = 1.5, y = 10.02, label = paste("p =", format(t_test_WT_KOWT$p.value, digits = 2)), size = 3) +
  annotate("text", x = 2.5, y = 10.025, label = paste("p =", format(t_test_KOWT_KONLS$p.value, digits = 2)), size = 3) +
  annotate("text", x = 3, y = 10.11, label = paste("p =", format(t_test_KOWT_KOSIG$p.value, digits = 2)), size = 3) +  
  annotate("text", x = 4, y = 10.21, label = paste("p =", format(t_test_KOWT_KONOG$p.value, digits = 2)), size = 3)


#nucleoside transporter
#Check ensemble ID for
#SLC29A1: ENSG00000112759
#SLC29A2: ENSG00000174669

gene_symbols <- c("ENSG00000112759", "ENSG00000174669")
filtered_counts <- rld_counts[rownames(count_data) %in% gene_symbols, ]
filtered_counts_2<- t(filtered_counts)
filtered_counts_2 <- as.data.frame(filtered_counts_2)
filtered_counts_2$condition <- rownames(filtered_counts_2)
filtered_counts_2 <- filtered_counts_2 %>%
  mutate(condition = sub("[0-9]+$", "", condition))
filtered_counts_2 <- filtered_counts_2[!filtered_counts_2$condition %in% c("K0NES"),]

SLC29A1 <- filtered_counts_2[,c(1,3)]
SLC29A2 <- filtered_counts_2[,c(2,3)]

ggplot(SLC29A1, aes(condition, ENSG00000112759, fill=condition )) + 
  geom_bar(stat="identity", width = 0.5) + 
  labs(x = "", y = "normalised gene count", title = "SLC29A1") + 
  theme_classic() +  
  scale_fill_grey() + theme(legend.position = "none")

ggplot(SLC29A2, aes(condition, ENSG00000174669, fill=condition )) + 
  geom_bar(stat="identity", width = 0.5) + 
  labs(x = "", y = "normalised gene count", title = "SLC29A2") + 
  theme_classic() +  scale_fill_grey() +theme(legend.position = "none")
