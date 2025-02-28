#################################
## Title: Extended Data Figure 6 panel e in Nomura et al - 3-axis per state
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: 3-axis per state
#################################

ggplot(feature_scores, aes(x = LineagePlot, y = C5)) +
  geom_point(aes(fill = State), color = "black", shape = 22, size = 3, alpha = .75) +
  geom_smooth(method = "loess", color = "black", se = F, size = 2, span = .75) +
  scale_fill_brewer(palette = "Set3") +
  xlab("ECM <-> Neuronal") +
  ylab("-> Glial") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())
