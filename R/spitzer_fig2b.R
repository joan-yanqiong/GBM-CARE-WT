#################################
## Title: Spitzer Figure 2 panel b - summary of assignment to composition, malignant states and baseline profile groups
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Produce barplots summarizing the assignment to composition, malignant states and baseline profile groups
#################################

####################################################################################################################################
# Composition panel
####################################################################################################################################

comp_cluster_stats_summarized <- comp_cluster_stats %>%
  group_by(SampleCluster, PvsR) %>%
  filter(!duplicated(Sample)) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PvsR) %>%
  mutate(N = sum(n), Freq = n / N, label = paste0(round(Freq * 100, 1), "%"), label_full = paste0(SampleCluster, " (", label, ")")) %>%
  ungroup()

ggplot(comp_cluster_stats_summarized, aes(x = SampleCluster, y = Freq, fill = PvsR)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_brewer(name = "", palette = "Set1") +
  # scale_y_continuous(labels = percent) +
  xlab("Composition group") +
  ylab("Proportion") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 16),
        panel.grid.major = element_line())

comp_cluster_stats_summarized %>%
  group_by(SampleCluster) %>%
  summarise(pval = fisher.test(x = matrix(data = c(c(n), c(N - n)), nrow = 2))$p.value)

comp_cluster_stats_summarized %>%
  filter(SampleCluster %in% c("LP - OC", "LP - GN")) %>%
  group_by(PvsR) %>%
  summarise(n = sum(n), N = first(N), .groups = "drop") %>%
  summarise(pval = fisher.test(x = matrix(data = c(c(n), c(N - n)), nrow = 2))$p.value)

####################################################################################################################################
# Malignant state panel
####################################################################################################################################

mal_cluster_stats_summarized <- mal_cluster_stats %>%
  group_by(SampleCluster, PvsR) %>%
  filter(!duplicated(Sample)) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PvsR) %>%
  mutate(N = sum(n), Freq = n / N, label = paste0(round(Freq * 100, 1), "%"), label_full = paste0(SampleCluster, " (", label, ")")) %>%
  ungroup()

ggplot(mal_cluster_stats_summarized, aes(x = SampleCluster, y = Freq, fill = PvsR)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent) +
  xlab("Malignant group") +
  ylab("Proportion") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 16),
        panel.grid.major = element_line())

mal_cluster_stats_summarized %>%
  group_by(SampleCluster) %>%
  summarise(pval = fisher.test(x = matrix(data = c(c(n), c(N - n)), nrow = 2))$p.value)

####################################################################################################################################
# Baseline profile panel
####################################################################################################################################

scp_stats_summarized <- dm %>%
  group_by(SCP, PvsR) %>%
  filter(!duplicated(ID)) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PvsR) %>%
  mutate(N = sum(n), Freq = n / N, label = paste0(round(Freq * 100, 1), "%"), label_full = paste0(SCP, " (", label, ")")) %>%
  ungroup()

ggplot(scp_stats_summarized, aes(x = SCP, y = Freq, fill = PvsR)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent) +
  xlab("BP group") +
  ylab("Proportion") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 16),
        panel.grid.major = element_line())

scp_stats_summarized %>%
  group_by(SCP) %>%
  summarise(pval = fisher.test(x = matrix(data = c(c(n), c(N - n)), nrow = 2))$p.value)
