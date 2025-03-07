#################################
## Title: Extended Data Figure 5d in Spitzer et al - cell type proportions across timepoints stratified by MGMT expression level
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot of cell type proportions across timepoints stratified by MGMT expression level
#################################

d <- comp_cluster_stats %>%
  filter(Patient %in% mgmt_md$Patient, Timepoint != "T3") %>%
  left_join(mgmt_md %>%
              select(Sample, MGMT, MGMT_exp_level),
            by = "Sample")

d_stats_summary <- acast(data = d, formula = Sample ~ CellType, value.var = "Freq")
d_stats_summary[is.na(d_stats_summary)] <- 0
d_stats_summary <- melt(d_stats_summary) %>%
  as_tibble()
colnames(d_stats_summary) <- c("Sample", "CellType", "Freq")

d_stats_summary <- d_stats_summary %>%
  left_join(d %>%
              ungroup() %>%
              select(Sample, Patient, Timepoint, PvsR, MGMT, MGMT_exp_level) %>%
              filter(!duplicated(Sample)),
            by = "Sample")

d_stats_summary <- d_stats_summary %>%
  group_by(Patient, CellType) %>%
  mutate(MGMT_exp_level = MGMT_exp_level[Timepoint == "T1"]) %>%
  ungroup()

pval_tbl <- d_stats_summary %>%
  group_by(CellType, MGMT_exp_level) %>%
  summarise(pval = wilcox.test(Freq[Timepoint == "T1"], Freq[Timepoint == "T2"])$p.value, .groups = "drop")

pval_tbl <- d_stats_summary %>%
  filter(PvsR == "Primary") %>%
  group_by(CellType) %>%
  summarise(pval = t.test(Freq[MGMT_exp_level == "High"],
                          Freq[MGMT_exp_level == "Low"], var.equal = T)$p.value, .groups = "drop") %>%
  mutate(padj = p.adjust(pval, "holm"))

d_stats_summary %>%
  filter(MGMT_exp_level %in% c("High", "Low"), CellType %ni% c("OPC", "Other")) %>%
  ggboxplot(x = "CellType", y = "Freq", add = "jitter", color = "MGMT_exp_level") +
  facet_grid(rows = vars(PvsR), scales = "free_y") +
  scale_color_brewer(name = "MGMT expression level", palette = "Set1") +
  scale_y_continuous(labels = percent) +
  xlab("Cell type") +
  ylab("Proportion") +
  theme_gbm_pvsr(legend.position = "top", axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())
