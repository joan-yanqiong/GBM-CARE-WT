#################################
## Title: Figure 3 panel c in Nomura et al - score difference between BP-Neuronal and BP-ECM
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: score difference between BP-Neuronal and BP-ECM to show the low intra-sample variability compared with inter-sample variability
#################################

feature_scores$Ord <- factor(feature_scores$ID,
                             feature_scores %>%
                               arrange(Diff) %>%
                               pull(ID))

gene_ord <- cor(t(pbm_list[pca_gene, ]), feature_scores$Diff) %>%
  as_tibble(rownames = "Gene") %>%
  arrange(V1)

sample_ord <- feature_scores %>%
  group_by(SID) %>%
  summarise(Mean = mean(Diff)) %>%
  arrange(Mean)

feature_scores$Sample_ord <- factor(feature_scores$SID, sample_ord$SID)

ft_mn <- feature_scores %>%
  group_by(Sample_ord) %>%
  summarise(Mean = mean(Diff), .groups = "drop")

ggplot(feature_scores, aes(x = Sample_ord, y = Diff)) +
  geom_point(data = ft_mn, mapping = aes(x = Sample_ord, y = Mean), shape = 21, color = "black", fill = "black", size = 3, alpha = .75, inherit.aes = F) +
  geom_line(data = ft_mn %>% mutate(grp = 1), mapping = aes(x = Sample_ord, y = Mean, group = grp), color = "black", size = 1, inherit.aes = F) +
  geom_line(data = feature_scores, mapping = aes(x = Sample_ord, y = Diff, group = Sample_ord), color = "black", size = .5, inherit.aes = F) +
  geom_point(aes(fill = State), color = "black", shape = 21, size = 2, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Samples") +
  ylab("Score difference (SCP-ECM - SCP-Neuronal)") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank())
