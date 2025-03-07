#################################
## Title: Spitzer Figure 4 panel b - MGMT expression in TME
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a boxplot of MGMT expression across TME cell types
#################################

# Compute MGMT expression level per sample
mgmt_exp_tme_tbl <- lapply(mgmt_md$Sample[mgmt_md$PvsR == "Primary"], function(sname) {
  print(sname)
  cts <- meta_data %>%
    filter(Sample == sname) %>%
    filter(CellType %in% c("Oligodendrocyte", "Macrophage", "Astrocyte", "Tcell",
                           "Excitatory neuron", "Inhibitory neuron", "Endothel", "Pericyte"))
  valid_cts <- cts %>%
    group_by(CellType) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n >= 20) %>%
    pull(CellType)
  
  res <- lapply(valid_cts, function(ct) {
    print(ct)
    cellids <- cts %>%
      filter(CellType == ct) %>%
      pull(CellID)
    exp <- umi_data_all[[sname]][, cellids]
    exp <- umi2upm(exp)
    tibble(Sample = sname, CellType = ct, Exp = log2(mean(exp["MGMT", ]) + 1))
  })
  res <- do.call(rbind, res)
  return(res)
})
mgmt_exp_tme_tbl <- do.call(rbind, mgmt_exp_tme_tbl)

mgmt_exp_tme_tbl <- mgmt_exp_tme_tbl %>%
  left_join(mgmt_md %>%
              select(Sample, MGMT), by = "Sample")

mgmt_exp_tme_tbl %>%
  filter(MGMT != "NOS") %>%
  group_by(CellType) %>%
  summarise(MeanMET = mean(Exp[MGMT == "MET"]), MeanUM = mean(Exp[MGMT == "UM"]),
            pval = wilcox.test(Exp[MGMT == "MET"], Exp[MGMT == "UM"])$p.value, .groups = "drop") %>%
  mutate(padj = p.adjust(pval, "holm"))

mgmt_exp_tme_mal_tbl <- rbind(mgmt_exp_tme_tbl,
                              mgmt_md %>%
                                filter(PvsR == "Primary") %>%
                                select(Sample, MGMT_exp, MGMT) %>%
                                mutate(CellType = "Malignant") %>%
                                rename(Exp = MGMT_exp) %>%
                                select(Sample, CellType, Exp, MGMT))

mgmt_exp_tme_mal_tbl %>%
  filter(MGMT != "NOS") %>%
  ggboxplot(x = "CellType", y = "Exp", color = "MGMT", add = "jitter",
            xlab = "Cell type", ylab = "MGMT expression [log2]") +
  scale_color_manual(name = "MGMT status", values = mgmt_color_vec) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())
