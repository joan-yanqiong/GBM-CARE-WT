#################################
## Title: Extended Data Figure 5h in Spitzer et al - differential gene expression across timepoints and MGMT expression levels
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a scatter plot showing the longitudinal gene expression difference across MGMT expression groups
#################################

####################################################################################################################################
# Generate the pseudobulk profiles and select genes for analysis
####################################################################################################################################

d <- mdata %>%
  filter(Sample %in% valid_samples$Sample)

mgmt_state_uncontrolled_m_all <- lapply(unique(d$Sample), function(sname) {
  print(sname)
  x <- umi_data_all[[sname]][, d$CellID[d$Sample == sname]]
  x <- umi2upm(x)
  x <- x[rownames(x) %ni% junk_genes, ]
  log(rowMeans(x) + 1)
})
mgmt_state_uncontrolled_m_all <- do.call(cbind, mgmt_state_uncontrolled_m_all)
colnames(mgmt_state_uncontrolled_m_all) <- unique(d$Sample)

rm <- rowMeans(mgmt_state_uncontrolled_m_all)

mgmt_state_uncontrolled_m_all <- mgmt_state_uncontrolled_m_all[names(rm[rm > 1]), ]

mgmt_state_uncontrolled_m <- mgmt_state_uncontrolled_m_all[, mgmt_md$Sample]

####################################################################################################################################
####################################################################################################################################
# T1 vs. T2 DEA per MGMT status
####################################################################################################################################
####################################################################################################################################

mgmt_t2_vs_t1 <- tibble(Gene = rownames(mgmt_state_uncontrolled_m),
                        MET = rowMeans(mgmt_state_uncontrolled_m[, mgmt_md$Sample[mgmt_md$Timepoint == "T2" & mgmt_md$MGMT_exp_level == "Low"]]) - rowMeans(mgmt_state_uncontrolled_m[, mgmt_md$Sample[mgmt_md$Timepoint == "T1" & mgmt_md$MGMT_exp_level == "Low"]]),
                        UM = rowMeans(mgmt_state_uncontrolled_m[, mgmt_md$Sample[mgmt_md$Timepoint == "T2" & mgmt_md$MGMT_exp_level == "High"]]) - rowMeans(mgmt_state_uncontrolled_m[, mgmt_md$Sample[mgmt_md$Timepoint == "T1" & mgmt_md$MGMT_exp_level == "High"]]))

mgmt_t2_vs_t1$Diff <- mgmt_t2_vs_t1$MET - mgmt_t2_vs_t1$UM

mgmt_t2_vs_t1 <- mgmt_t2_vs_t1 %>%
  group_by(Gene) %>%
  mutate(Dist = abs(UM - MET) / sqrt(2) * sign(MET - UM)) %>%
  ungroup()

min_diff_th <- log2(1.5)

met_high_genes <- mgmt_t2_vs_t1 %>%
  filter(MET > 0 & UM < 0) %>%
  arrange(desc(Dist)) %>%
  head(100) %>%
  pull(Gene)

um_high_genes <- mgmt_t2_vs_t1 %>%
  filter(UM > 0 & MET < 0) %>%
  arrange(Dist) %>%
  head(100) %>%
  pull(Gene)
length(met_high_genes)
length(um_high_genes)

rec_up_cons <- mgmt_t2_vs_t1 %>%
  filter(MET > .25*min_diff_th & UM > .25*min_diff_th & Gene %ni% unique(c(met_high_genes, um_high_genes))) %>%
  arrange(desc(MET + UM)) %>%
  head(100) %>%
  pull(Gene)

rec_dn_cons <- mgmt_t2_vs_t1 %>%
  filter(MET < -.25*min_diff_th & UM < -.25*min_diff_th & Gene %ni% unique(c(met_high_genes, um_high_genes))) %>%
  arrange(MET + UM) %>%
  head(100) %>%
  pull(Gene)
length(rec_up_cons)
length(rec_dn_cons)

mgmt_t2_vs_t1 <- mgmt_t2_vs_t1 %>%
  mutate(GeneGroup = case_when(Gene %in% met_high_genes ~ "MGMT_Low",
                               Gene %in% um_high_genes ~ "MGMT_High",
                               Gene %in% rec_up_cons ~ "CONS_UP",
                               Gene %in% rec_dn_cons ~ "CONS_DN",
                               TRUE ~ "NA")) %>%
  mutate(GeneGroup = factor(GeneGroup, c("MGMT_Low", "MGMT_High", "CONS_UP", "CONS_DN", "NA")))
table(mgmt_t2_vs_t1$GeneGroup)

####################################################################################################################################
# Plot panel
####################################################################################################################################

mgmt_sig_color_vec <- c("MGMT_Low" = "orange", "MGMT_High" = "darkgreen", "CONS_UP" = "purple", "CONS_DN" = "dodgerblue", "NA" = "grey")

ggscatter(data = mgmt_t2_vs_t1, x = "UM", y = "MET", xlab = "MGMT-High T2 vs. T1", ylab = "MGMT-Low T2 vs. T1", color = "GeneGroup") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 2) +
  scale_color_manual(name = "Gene group", values = mgmt_sig_color_vec) +
  scale_x_continuous(limits = c(-1.5, 1.5), oob = squish, breaks = seq(-1.5, 1.5, .5)) +
  scale_y_continuous(limits = c(-1.5, 1.5), oob = squish, breaks = seq(-1.5, 1.5, .5)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())
