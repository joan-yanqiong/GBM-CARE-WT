#################################
## Title: Extended Data Figure 2a in Spitzer et al - composition differences across timepoints
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot of composition differences across timepoints
#################################

d_stats_summary <- meta_data %>%
  filter(Timepoint != "T3") %>%
  group_by(Sample, CellType_reduced) %>%
  summarize(n = n()) %>%
  group_by(Sample) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup()

d_stats_summary <- acast(data = d_stats_summary, formula = Sample ~ CellType_reduced, value.var = "Freq")
d_stats_summary[is.na(d_stats_summary)] <- 0
d_stats_summary <- melt(d_stats_summary) %>%
  as_tibble()
colnames(d_stats_summary) <- c("Sample", "CellType", "Freq")
dim(d_stats_summary)

d_stats_summary <- d_stats_summary %>%
  left_join(sample_data %>%
              ungroup() %>%
              dplyr::select(Sample, Patient, Timepoint, PvsR),
            by = "Sample")

d_stats_summary <- d_stats_summary %>%
  mutate(PvsR = as.character(PvsR)) %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrence"))

d_stats_summary %>%
  ggplot(aes(x = Timepoint, y = Freq, fill = CellType)) +
  facet_grid(cols = vars(CellType), scales = "free") +
  geom_boxplot() +
  # geom_jitter(aes(group = Timepoint), shape = 21, color = "black", fill = "black", alpha = .5, width = .2) +
  # geom_point(aes(group = Patient), position = position_jitter(seed = 1), shape = 21, color = "black", fill = "black", alpha = .5) +
  # geom_line(aes(group = Timepoint), size = 1, linetype = "dashed", color = "black", alpha = .5) +
  geom_line(aes(group = Patient), size = .5, linetype = "dashed", color = "black", alpha = .5) +
  scale_fill_brewer(name = "", palette = "Spectral") +
  xlab("") +
  ylab("Proportion") +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "none", axis.text.x = element_text(size = 12)) +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))

####################################################################################################################################
# Statistics
####################################################################################################################################

d_stats_summary <- meta_data %>%
  group_by(Sample, CellType_reduced) %>%
  summarize(n = n()) %>%
  group_by(Sample) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup()

d_stats_summary <- acast(data = d_stats_summary, formula = Sample ~ CellType_reduced, value.var = "Freq")
d_stats_summary[is.na(d_stats_summary)] <- 0
d_stats_summary <- melt(d_stats_summary) %>%
  as_tibble()
colnames(d_stats_summary) <- c("Sample", "CellType", "Freq")

d_stats_summary <- d_stats_summary %>%
  left_join(sample_data %>%
              ungroup() %>%
              select(Sample, Patient, Timepoint, PvsR),
            by = "Sample")

d_stats_summary <- d_stats_summary %>%
  mutate(PvsR = as.character(PvsR)) %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrence"))

d_stats_summary %>%
  filter(CellType == "Malignant") %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

wilcox.test(x = d_stats_summary %>% filter(CellType == "Malignant", PvsR == "Primary") %>% pull(Freq),
            y = d_stats_summary %>% filter(CellType == "Malignant", PvsR == "Recurrence") %>% pull(Freq), alternative = "greater")

d_stats_summary %>%
  filter(CellType == "Oligodendrocyte") %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

wilcox.test(x = d_stats_summary %>% filter(CellType == "Oligodendrocyte", PvsR == "Recurrence") %>% pull(Freq),
            y = d_stats_summary %>% filter(CellType == "Oligodendrocyte", PvsR == "Primary") %>% pull(Freq), alternative = "greater")

d_stats_summary %>%
  filter(CellType %in% c("Excitatory neuron", "Inhibitory neuron")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

d_stats_summary %>%
  filter(CellType %in% c("Excitatory neuron", "Inhibitory neuron")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  ungroup() %>%
  summarise(pval = wilcox.test(x = Freq[PvsR == "Recurrence"], y = Freq[PvsR == "Primary"], alternative = "greater")$p.value)

d_stats_summary %>%
  filter(CellType == "Astrocyte") %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

wilcox.test(x = d_stats_summary %>% filter(CellType == "Astrocyte", PvsR == "Recurrence") %>% pull(Freq),
            y = d_stats_summary %>% filter(CellType == "Astrocyte", PvsR == "Primary") %>% pull(Freq), alternative = "greater")
