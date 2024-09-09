title: "Cell_death"
author: "Alisa Schmidt"

#plot data with all values for recon, paired t test
data <- read.delim("/Users/alisc/Desktop/CRG/Cell_death/all_recon.csv", sep=",")
data <- data[,1:5]


data_mean <- data %>%
  group_by(cell, drug, death.type) %>%
  mutate(mean = mean(value1),
         sd_1 = sd(value1)) %>%
  ungroup()


t_test_results <- data %>%
  pivot_wider(names_from = drug, values_from = value1) %>%
  group_by(death.type, cell) %>%
  summarise(p_value = t.test(DMSO, Etoposide, paired = TRUE)$p.value) %>%
  ungroup() %>%
  mutate(p.signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  ))



ordered <- c("WT-empty","KO-empty","KO-wt","KO-nls")
data_mean$cell <- factor(data$cell, levels = ordered)

ggplot(data_mean, aes(x = cell, y = value1, color = drug)) +
  geom_point(position = position_dodge(width = 0.8)) +  # Points for value1
  facet_wrap(~death.type) +  # Separate plots based on death_type
  geom_bar(aes(y = mean, group=drug), stat = "identity", position = position_dodge(width = 0.9), color = "black", alpha = 0) +  # Bars for mean, black color
  geom_errorbar(aes(ymin = mean - sd_1, ymax = mean + sd_1, group = interaction(drug, death.type)), 
                position = position_dodge(width = 0.8), width = 0.25, color = "black") + 
  theme_minimal() +
  scale_color_manual(values = c("DMSO" = "blue", "Etoposide" = "red")) +  # Custom colors for points
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
