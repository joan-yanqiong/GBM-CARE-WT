#################################
## Title: Spitzer Figure 4 panel g - signatures associated with MGMT high and low expression levels
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a heatmap showing the signatures associated with MGMT high and low expression levels
#################################

####################################################################################################################################
# Generate a DE vector for each patient
####################################################################################################################################

mgmt_t2_vs_t1_state_uncontrolled_m_all <- lapply(unique(mgmt_md$Patient), function(pt) {
  t1 <- mgmt_md$Sample[mgmt_md$Patient == pt & mgmt_md$Timepoint == "T1"]
  t2 <- mgmt_md$Sample[mgmt_md$Patient == pt & mgmt_md$Timepoint == "T2"]
  mgmt_state_uncontrolled_m_all[, t2] - mgmt_state_uncontrolled_m_all[, t1]
})
mgmt_t2_vs_t1_state_uncontrolled_m_all <- do.call(cbind, mgmt_t2_vs_t1_state_uncontrolled_m_all)
dim(mgmt_t2_vs_t1_state_uncontrolled_m_all)
colnames(mgmt_t2_vs_t1_state_uncontrolled_m_all) <- unique(mgmt_md$Patient)

mgmt_t2_vs_t1_state_uncontrolled_m <- mgmt_t2_vs_t1_state_uncontrolled_m_all[, mgmt_md$Patient]

####################################################################################################################################
# Generate per-patient meta-data to hold clinical data and signature scores
####################################################################################################################################

mgmt_md_p_all <- mgmt_md %>%
  group_by(Patient) %>%
  summarise(SurgicalInterval = first(SurgicalInterval), SurgicalInterval2 = first(SurgicalInterval2),
            OSmonth = first(OSmonth), VT = first(VT), Age = first(Age), Gender = first(Gender),
            SCP_T1 = SCP[Timepoint == "T1"], SCP = SCP[Timepoint == "T2"], #SCP2 = SCP2[Timepoint == "T2"],
            MGMT = MGMT[Timepoint == "T1"],
            MGMT_exp = MGMT_exp[Timepoint == "T2"] - MGMT_exp[Timepoint == "T1"],
            MGMT_exp_level = MGMT_exp_level[Timepoint == "T1"],
            # dsDNA_break_repair_score = dsDNA_break_repair_score[Timepoint == "T2"] - dsDNA_break_repair_score[Timepoint == "T1"],
            .groups = "drop")

mgmt_md_p_all$MGMT_MET <- colMeans(mgmt_t2_vs_t1_state_uncontrolled_m_all[met_high_genes, mgmt_md_p_all$Patient])
mgmt_md_p_all$MGMT_UM <- colMeans(mgmt_t2_vs_t1_state_uncontrolled_m_all[um_high_genes, mgmt_md_p_all$Patient])
mgmt_md_p_all$REC_UP <- colMeans(mgmt_t2_vs_t1_state_uncontrolled_m_all[rec_up_cons, mgmt_md_p_all$Patient])
mgmt_md_p_all$REC_DN <- colMeans(mgmt_t2_vs_t1_state_uncontrolled_m_all[rec_dn_cons, mgmt_md_p_all$Patient])

mgmt_md_p_all$MGMT_Diff <- mgmt_md_p_all$MGMT_MET - mgmt_md_p_all$MGMT_UM
mgmt_md_p_all$CONS_Diff <- mgmt_md_p_all$REC_UP - mgmt_md_p_all$REC_DN

mgmt_md_p <- mgmt_md_p_all %>%
  filter(MGMT_exp_level != "Intermediate")

####################################################################################################################################
# Plot heatmap of MET vs. UM sigs
####################################################################################################################################

d <- mgmt_md_p

dm <- mgmt_t2_vs_t1_state_uncontrolled_m[c(met_high_genes, um_high_genes), d$Patient]
dim(dm)

gene_ord <- cor(t(dm), d$MGMT_Diff)
gene_ord <- setNames(gene_ord[, 1], rownames(gene_ord))
gene_ord <- sort(gene_ord)

pt_ord <- d %>%
  arrange(MGMT_Diff)

dm <- melt(dm) %>%
  as_tibble()

dm$Var1 <- factor(as.character(dm$Var1), names(gene_ord))
dm$Var2 <- factor(as.character(dm$Var2), pt_ord$Patient)

ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(name = "log2FC", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-1, 1), oob = squish, labels = c("-2", "", "0", "", "2")) +
  xlab("Genes") +
  ylab("Patients") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 # axis.text.x = element_blank(),
                 axis.text.x = element_text(size = 5, angle = 90),
                 # axis.ticks.x = element_blank(),
                 axis.text.y = element_text(size = 12))

gene_sig_dm <- tibble(Gene = names(gene_ord)) %>%
  mutate(Sig = case_when(Gene %in% met_high_genes ~ "MGMT_Low",
                         Gene %in% um_high_genes ~ "MGMT_High")) %>%
  mutate(Gene = factor(Gene, Gene))

ggplot(gene_sig_dm, aes(x = Gene, y = "Sig", fill = Sig)) +
  geom_tile() +
  scale_fill_manual(values = mgmt_sig_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_blank())

mgmt_ord_dm <- tibble(Patient = pt_ord$Patient,
                      MGMT = pt_ord$MGMT_group) %>%
  mutate(Patient = factor(Patient, Patient))

mgmt_exp_color_vec <- c("High" = "red", "Low" = "black", "Intermediate" = "grey")

ggplot(mgmt_ord_dm, aes(x = Patient, y = "MGMT", fill = MGMT)) +
  geom_tile() +
  scale_fill_manual(values = mgmt_exp_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 10))
