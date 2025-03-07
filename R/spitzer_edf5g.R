#################################
## Title: Extended Data Figure 5g in Spitzer et al - baseline profile difference across timepoints stratified by MGMT expression level
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot of baseline profile difference across timepoints stratified by MGMT expression level
#################################

d <- feature_scores_all_per_sample %>%
  select(Sample, Patient, Timepoint, PvsR, C1, C3, C5, SCP) %>%
  filter(Patient %in% mgmt_md$Patient) %>%
  left_join(mgmt_md %>%
              select(Sample, MGMT, MGMT_exp_level),
            by = "Sample")

dm <- melt(d) %>%
  as_tibble()

d_stats_summary_t2_vs_t1 <- dm %>%
  group_by(variable, Patient, MGMT_exp_level) %>%
  summarise(Diff = value[Timepoint == "T2"] - value[Timepoint == "T1"], .groups = "drop")

pval_tbl <- d_stats_summary_t2_vs_t1 %>%
  group_by(variable) %>%
  summarise(MeanLow = mean(Diff[MGMT_exp_level == "Low"]),
            MeanHigh = mean(Diff[MGMT_exp_level == "High"]),
            pval = t.test(Diff[MGMT_exp_level == "High"], Diff[MGMT_exp_level == "Low"], var.equal = T)$p.value,
            .groups = "drop") %>%
  mutate(padj = p.adjust(pval, "fdr"))

d_stats_summary_t2_vs_t1 %>%
  filter(MGMT_exp_level != "Intermediate") %>%
  ggboxplot(x = "variable", y = "Diff", add = "jitter", color = "MGMT_exp_level") +
  scale_color_brewer(name = "MGMT status", palette = "Set1") +
  scale_x_discrete(labels = c("C1" = "BP-ECM", "C3" = "BP-Neuronal", "C5" = "BP-Glial")) +
  xlab("BP") +
  ylab("Score") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top",
                 axis.text.x = element_text(size = 16)) +
  theme(panel.grid.major = element_line())
