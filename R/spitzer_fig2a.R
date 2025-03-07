#################################
## Title: Spitzer Figure 2 panel a - significant composition changes
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a boxplot of significant composition changes across timepoints
#################################

d_stats_summary <- d_stats_summary %>%
  filter(CellType %in% c("Malignant", "Oligodendrocyte", "Astrocyte", "Excitatory neuron", "Inhibitory neuron")) %>%
  mutate(CellType = as.character(CellType),
         CellType = ifelse(CellType %in% c("Excitatory neuron", "Inhibitory neuron"), "Neuron", CellType),
         CellType = factor(CellType, c("Malignant", "Oligodendrocyte", "Astrocyte", "Neuron")))

ggboxplot(d_stats_summary, x = "CellType", y = "Freq", add = "jitter", color = "PvsR") +
  scale_color_brewer(name = "", palette = "Set2") +
  scale_y_continuous(labels = percent) +
  xlab("Cell type") +
  ylab("Proportion") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top", axis.text.x = element_text(size = 16, angle = 0)) +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))
