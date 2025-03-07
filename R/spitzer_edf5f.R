#################################
## Title: Extended Data Figure 5f in Spitzer et al - cell type proportion difference across timepoints stratified by MGMT expression level
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot of cell type proportion difference across timepoints stratified by MGMT expression level
#################################

d_stats_summary_t2_vs_t1 <- d_stats_summary %>%
  group_by(CellType, Patient, MGMT_exp_level) %>%
  summarise(Diff = Freq[Timepoint == "T2"] - Freq[Timepoint == "T1"], .groups = "drop")

pval_tbl <- d_stats_summary_t2_vs_t1 %>%
  group_by(CellType) %>%
  summarise(MeanLow = median(Diff[MGMT_exp_level == "Low"]),
            MeanHigh = median(Diff[MGMT_exp_level == "High"]),
            pval = wilcox.test(Diff[MGMT_exp_level == "High"], Diff[MGMT_exp_level == "Low"])$p.value,
            .groups = "drop") %>%
  mutate(padj = p.adjust(pval, "fdr"))

pval_tbl <- d_stats_summary_t2_vs_t1 %>%
  group_by(CellType) %>%
  summarise(MeanLow = mean(Diff[MGMT_exp_level == "Low"]),
            MeanHigh = mean(Diff[MGMT_exp_level == "High"]),
            pval = t.test(Diff[MGMT_exp_level == "High"], Diff[MGMT_exp_level == "Low"], var.equal = T)$p.value,
            .groups = "drop") %>%
  mutate(padj = p.adjust(pval, "fdr"))

d_stats_summary_t2_vs_t1 %>%
  filter(MGMT_exp_level != "Intermediate", CellType %ni% c("OPC", "Other")) %>%
  ggboxplot(x = "CellType", y = "Diff", add = "jitter", color = "MGMT_exp_level") +
  scale_color_brewer(name = "MGMT status", palette = "Set1") +
  scale_y_continuous(labels = percent, limits = c(-.2, .2), oob = squish) +
  xlab("Cell type") +
  ylab("Proportion") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top",
                 axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())
