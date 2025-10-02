#################################
## Title: Extended Data Figure 3a in Spitzer et al - malignant state differences across timepoints
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot of malignant state differences across timepoints
#################################

valid_samples <- mdata %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n >= 50) %>%
  filter(Timepoint != "T3") %>%
  group_by(Patient) %>%
  filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
  ungroup()
dim(valid_samples)
table(valid_samples$Patient) %>% sort()

d_stats_summary <- mdata %>%
  filter(Sample %in% valid_samples$Sample) %>%
  group_by(Sample, State) %>%
  summarize(n = n()) %>%
  group_by(Sample) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup()

d_stats_summary <- acast(data = d_stats_summary, formula = Sample ~ State, value.var = "Freq")
d_stats_summary[is.na(d_stats_summary)] <- 0
d_stats_summary <- melt(d_stats_summary) %>%
  as_tibble()
colnames(d_stats_summary) <- c("Sample", "State", "Freq")

d_stats_summary <- d_stats_summary %>%
  left_join(sample_data %>%
              ungroup() %>%
              dplyr::select(Sample, Patient, Timepoint, PvsR),
            by = "Sample")

d_stats_summary <- d_stats_summary %>%
  mutate(PvsR = as.character(PvsR)) %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrence")) %>%
  mutate(PvsR = factor(PvsR, c("Primary", "Recurrence")))

d_stats_summary %>%
  mutate(State = factor(as.character(State), c("Cilia", "AC", "MES", "Hypoxia", "Stress", "GPC", "OPC", "NPC", "Neuron", "Cycling", "Unresolved"))) %>%
  mutate(Pair = Timepoint) %>%
  ggplot(aes(x = Timepoint, y = Freq, fill = State)) +
  facet_grid(cols = vars(State), scales = "free") +
  geom_boxplot() +
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

d_stats_summary <- mdata %>%
  filter(Sample %in% valid_samples$Sample) %>%
  group_by(Sample, State) %>%
  summarize(n = n()) %>%
  group_by(Sample) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup()

d_stats_summary <- acast(data = d_stats_summary, formula = Sample ~ State, value.var = "Freq")
d_stats_summary[is.na(d_stats_summary)] <- 0
d_stats_summary <- melt(d_stats_summary) %>%
  as_tibble()
colnames(d_stats_summary) <- c("Sample", "State", "Freq")

d_stats_summary <- d_stats_summary %>%
  left_join(sample_data %>%
              ungroup() %>%
              dplyr::select(Sample, Patient, Timepoint, PvsR),
            by = "Sample")

d_stats_summary <- d_stats_summary %>%
  mutate(PvsR = as.character(PvsR)) %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrence"))

pval_tbl <- d_stats_summary %>%
  group_by(State) %>%
  summarise(P = mean(Freq[Timepoint == "T1"]), R = mean(Freq[Timepoint == "T2"]),
            pval = wilcox.test(x = Freq[Timepoint == "T1"], y = Freq[Timepoint == "T2"], paired = T)$p.value, .groups = "drop")

d_stats_summary %>%
  filter(State %in% c("Neuron", "OPC")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

d_stats_summary %>%
  filter(State %in% c("Neuron", "OPC")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  ungroup() %>%
  summarise(pval = wilcox.test(x = Freq[PvsR == "Recurrence"], y = Freq[PvsR == "Primary"], alternative = "greater")$p.value)

d_stats_summary %>%
  filter(State %in% c("Neuron")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

d_stats_summary %>%
  filter(State %in% c("Neuron")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  ungroup() %>%
  summarise(pval = t.test(x = Freq[PvsR == "Recurrence"], y = Freq[PvsR == "Primary"], alternative = "greater")$p.value)

d_stats_summary %>%
  filter(State %in% c("AC", "MES", "GPC")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  group_by(PvsR) %>%
  summarise(Mean = mean(Freq))

d_stats_summary %>%
  filter(State %in% c("AC", "MES", "GPC")) %>%
  group_by(Sample, PvsR) %>%
  summarise(Freq = sum(Freq)) %>%
  ungroup() %>%
  summarise(pval = wilcox.test(x = Freq[PvsR == "Primary"], y = Freq[PvsR == "Recurrence"], alternative = "greater")$p.value)

pval_tbl <- d_stats_summary %>%
  group_by(State) %>%
  summarise(pval = wilcox.test(Freq[PvsR == "Primary"], Freq[PvsR != "PRimary"])$p.value, .groups = "drop")
