#################################
## Title: Extended Data Figure 7 panel b in Nomura et al - Correlation between snRNAseq and Spatial BPs
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a heatmap showing the Correlation between snRNAseq and Spatial BPs
#################################

bp_cor_tbl <- read.csv(file = paste0(TABLES_ROOT, "snRNAseq_Spatial_BP_cor.csv"),
                       header = T, stringsAsFactors = F) %>%
  as_tibble()
bp_cor_tbl <- bp_cor_tbl[-1, -5]
bp_cor_tbl$padj <- p.adjust(bp_cor_tbl$pval, "holm")

bp_cor_tbl$Spatial_BP <- factor(bp_cor_tbl$Spatial_BP, c("BP1", "BP3", "BP2"))

ggplot(bp_cor_tbl, aes(x = snRNAseq_BP, y = Spatial_BP, fill = cor)) +
  geom_tile(color = "black") +
  geom_text(aes(label = padj)) +
  scale_fill_gradient2(name = "Correlation", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-.5, .5), oob = squish) +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())
