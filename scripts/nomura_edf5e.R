#################################
## Title: Extended Data Figure 5 panel e in Nomura et al - observed vs. expected proportion of hybrids
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a barplot showing the observed vs. expected proportion of hybrids per hybrid pair
#################################

obs_vs_exp_hybrids$Label <- obs_vs_exp_hybrids$padj

ggplot(obs_vs_exp_hybrids, aes(x = Hybrid, y = log2FC_2, fill = sigRes)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = Label, y = 3), color = "black", size = 10) +
  scale_fill_brewer(name = "", palette = "Set1") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 1) +
  xlab("Hybrid pair") +
  ylab("Observed vs. expected\ntechnical hgybrids [log2]") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90), legend.position = "none") +
  theme(panel.grid.major = element_line())
