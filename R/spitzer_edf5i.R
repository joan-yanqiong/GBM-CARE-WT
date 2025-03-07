#################################
## Title: Extended Data Figure 5i in Spitzer et al - signatures up and downregulated regardless of MGMT expression levels
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Generate a heatmap showing the signatures up and downregulated regardless of MGMT expression levels
#################################

dm <- mgmt_t2_vs_t1_state_uncontrolled_m[c(rec_up_cons, rec_dn_cons), mgmt_md_p$Patient]
dim(dm)

mgmt_md_p$REC_Diff <- mgmt_md_p$REC_UP - mgmt_md_p$REC_DN

gene_ord <- cor(t(dm), mgmt_md_p$REC_Diff)
gene_ord <- setNames(gene_ord[, 1], rownames(gene_ord))
gene_ord <- sort(gene_ord)

pt_ord <- mgmt_md_p %>%
  arrange(REC_Diff)

dm <- melt(dm) %>%
  as_tibble()

dm$Var1 <- factor(as.character(dm$Var1), names(gene_ord))
dm$Var2 <- factor(as.character(dm$Var2), pt_ord$Patient)

ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(name = "log2FC", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-1, 1), oob = squish, labels = c("-2", "", "0", "", "2")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_text(size = 5, angle = 90),
                 # axis.text.x = element_blank(),
                 axis.text.y = element_text(size = 12))

gene_sig_dm <- tibble(Gene = names(gene_ord)) %>%
  mutate(Sig = case_when(Gene %in% rec_up_cons ~ "CONS_UP",
                         Gene %in% rec_dn_cons ~ "CONS_DN")) %>%
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

ggplot(mgmt_ord_dm, aes(x = Patient, y = "MGMT", fill = MGMT)) +
  geom_tile() +
  scale_fill_manual(values = mgmt_exp_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 10))
