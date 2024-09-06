title: "Gene_counts"
author: "Alisa Schmidt"
date: "07.06.2024"

library(DESeq2)
library(ggplot2)
library(tidyr)
library(tidyverse)

#load the DeSeq data 
load("/Users/alisc/Desktop/CRG/DeSeq R/DDS_file.RData")
count_data <- counts(dds)
rld <- rlog(dds, blind = F)
rld_counts <- assay(rld)

IMPDH2_counts <- rld_counts[rownames(count_data) == "ENSG00000178035", ]
IMPDH2_counts

NAMPT_counts <- rld_counts[rownames(count_data) == "ENSG00000105835", ]
NAMPT_counts

NAMPT_counts2 <- as.data.frame(NAMPT_counts)
NAMPT_counts2$condition <- rownames(NAMPT_counts2)


NAMPT_counts2 <- NAMPT_counts2 %>%
  mutate(condition = sub("[0-9]+$", "", condition))
NAMPT_counts2 <- NAMPT_counts2[!NAMPT_counts2$condition %in% c("K0NES", ("K0SIG"), ("K0N0G")),]
NAMPT_counts2_mean <- NAMPT_counts2 %>%
  group_by(condition) %>%
  summarise(mean(NAMPT_counts), sd(NAMPT_counts))
colnames(NAMPT_counts2_mean) <- c("condition", "Mean", "Std")

#t test
wt_data <- NAMPT_counts2 %>% filter(condition == "WT") %>% select(NAMPT_counts)
kowt_data <- NAMPT_counts2 %>% filter(condition == "K0WT") %>% select(NAMPT_counts)
konls_data <- NAMPT_counts2 %>% filter(condition == "K0NLS") %>% select(NAMPT_counts)
t_test_result <- t.test(wt_data$NAMPT_counts, kowt_data$NAMPT_counts, paired = FALSE)
t_test_result2 <- t.test(konls_data$NAMPT_counts, kowt_data$NAMPT_counts, paired = FALSE)

print(t_test_result)

order <- c("WT", "K0WT", "K0NLS")
NAMPT_counts2_mean$condition <- factor(NAMPT_counts2_mean$condition, levels = order)
#pdf("/Users/alisc/Desktop/CRG/gene_counts/NAMPT_counts.pdf")
ggplot(NAMPT_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean NAMPT Counts with Standard Deviation",
       y = "Mean NAMPT Counts")   + annotate("text", x = 1, y = max(NAMPT_counts2_mean$Mean) + 0.1,
                                             label = paste("p=", format(t_test_result$p.value, digits = 2)),
                                             size = 4, hjust = 0)  + annotate("text", x = 2, y = max(NAMPT_counts2_mean$Mean) + 0.1,
                                                                              label = paste("p=", format(t_test_result2$p.value, digits = 2)),
                                                                              size = 4, hjust = 0)
#dev.off()
my_comp <-  list(c("WT", "K0WT"), c("K0WT", "K0NLS"))

#pdf("/Users/alisc/Desktop/CRG/gene_counts/NAMPT_counts_t_test.pdf")
ggplot(NAMPT_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean NAMPT Counts with Standard Deviation",
       y = "Mean NAMPT Counts") + geom_segment(aes(x = 1, xend = 2, y = 9.49, yend = 9.49)) + # Line for WT vs KOWT
  geom_segment(aes(x = 2, xend = 3, y = 9.45, yend = 9.45)) + # Line for KOWT vs KONLS
  annotate("text", x = 1.5, y = 9.5, label = paste("p =", format(t_test_result$p.value, digits = 2)), size = 3) + # P-value for WT vs KOWT
  annotate("text", x = 2.5, y = 9.46, label = paste("p =", format(t_test_result2$p.value, digits = 2)), size = 3) 
#dev.off()

IMPDH2_counts2 <- as.data.frame(IMPDH2_counts)
IMPDH2_counts2$condition <- rownames(IMPDH2_counts2)


IMPDH2_counts2 <- IMPDH2_counts2 %>%
  mutate(condition = sub("[0-9]+$", "", condition))
IMPDH2_counts2 <- IMPDH2_counts2[!IMPDH2_counts2$condition %in% c("K0NES"),]
IMPDH2_counts2_mean <- IMPDH2_counts2 %>%
  group_by(condition) %>%
  summarise(mean(IMPDH2_counts), sd(IMPDH2_counts))
colnames(IMPDH2_counts2_mean) <- c("condition", "Mean", "Std")

ggplot(IMPDH2_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean IMPDH2 Counts with Standard Deviation",
       y = "Mean normalised IMPDH2 Counts")

IMPDH2_counts2_mean_recon <- IMPDH2_counts2_mean[!IMPDH2_counts2_mean$condition %in% c("K0NES", "K0N0G", "K0SIG"),]
IMPDH2_counts2_mean_recon$condition <- factor(IMPDH2_counts2_mean_recon$condition, levels = order)


ggplot(NAMPT_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean NAMPT Counts with Standard Deviation",
       y = "Mean NAMPT Counts")

ggplot(IMPDH2_counts2_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean IMPDH2 Counts with Standard Deviation",
       y = "Mean normalised IMPDH2 Counts")
