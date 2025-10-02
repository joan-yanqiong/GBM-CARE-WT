#################################
## Title: Extended Data Figure 6 panel c in Nomura et al - baseline profile scores across samples
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: plot scores of baseline profiles to show the low intra-sample variability compared with inter-sample variability
#################################

tmp <- feature_scores %>%
  mutate(Sample_ord = factor(Sample, feature_scores %>%
                               group_by(Sample) %>%
                               summarise(Mean = mean(C1)) %>%
                               arrange(Mean) %>%
                               pull(Sample)))
ft_mn <- tmp %>%
  group_by(Sample_ord) %>%
  summarise(Mean = mean(C1), .groups = "drop")

p1 <- ggplot(tmp, aes(x = Sample_ord, y = C1)) +
  geom_point(data = ft_mn, mapping = aes(x = Sample_ord, y = Mean), shape = 21, color = "black", fill = "black", size = 3, alpha = .75, inherit.aes = F) +
  geom_line(data = ft_mn %>% mutate(grp = 1), mapping = aes(x = Sample_ord, y = Mean, group = grp), color = "black", size = 1, inherit.aes = F) +
  geom_line(data = tmp, mapping = aes(x = Sample_ord, y = C1, group = Sample_ord), color = "black", size = .5, inherit.aes = F) +
  geom_point(aes(fill = State), color = "black", shape = 21, size = 2, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("BP-ECM score") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank())

tmp <- feature_scores %>%
  mutate(Sample_ord = factor(Sample, feature_scores %>%
                               group_by(Sample) %>%
                               summarise(Mean = mean(C3)) %>%
                               arrange(Mean) %>%
                               pull(Sample)))
ft_mn <- tmp %>%
  group_by(Sample_ord) %>%
  summarise(Mean = mean(C3), .groups = "drop")

p2 <- ggplot(tmp, aes(x = Sample_ord, y = C3)) +
  geom_point(data = ft_mn, mapping = aes(x = Sample_ord, y = Mean), shape = 21, color = "black", fill = "black", size = 3, alpha = .75, inherit.aes = F) +
  geom_line(data = ft_mn %>% mutate(grp = 1), mapping = aes(x = Sample_ord, y = Mean, group = grp), color = "black", size = 1, inherit.aes = F) +
  geom_line(data = tmp, mapping = aes(x = Sample_ord, y = C3, group = Sample_ord), color = "black", size = .5, inherit.aes = F) +
  geom_point(aes(fill = State), color = "black", shape = 21, size = 2, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("BP-Neuronal score") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank())

tmp <- feature_scores %>%
  mutate(Sample_ord = factor(Sample, feature_scores %>%
                               group_by(Sample) %>%
                               summarise(Mean = mean(C5)) %>%
                               arrange(Mean) %>%
                               pull(Sample)))
ft_mn <- tmp %>%
  group_by(Sample_ord) %>%
  summarise(Mean = mean(C5), .groups = "drop")

p3 <- ggplot(tmp, aes(x = Sample_ord, y = C5)) +
  geom_point(data = ft_mn, mapping = aes(x = Sample_ord, y = Mean), shape = 21, color = "black", fill = "black", size = 3, alpha = .75, inherit.aes = F) +
  geom_line(data = ft_mn %>% mutate(grp = 1), mapping = aes(x = Sample_ord, y = Mean, group = grp), color = "black", size = 1, inherit.aes = F) +
  geom_line(data = tmp, mapping = aes(x = Sample_ord, y = C5, group = Sample_ord), color = "black", size = .5, inherit.aes = F) +
  geom_point(aes(fill = State), color = "black", shape = 21, size = 2, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("BP-Glial score") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p1 + p2 + p3 + plot_layout(ncol = 1)
