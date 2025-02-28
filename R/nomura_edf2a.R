#################################
## Title: Extended Data Figure 2 panel a in Nomura et al - umap colored by generating lab
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Cohort-level umap colored by generating lab
#################################

d <- as_tibble(sds@meta.data)

umap_res <- Reductions(sds, "umap_corr")
d$UMAP1 <- umap_res@cell.embeddings[d$CellID, 1]
d$UMAP2 <- umap_res@cell.embeddings[d$CellID, 2]

ggplot(d, aes(x = UMAP1, y = UMAP2, color = Lab)) +
  ggrastr::geom_point_rast(size = .25) +
  scale_color_discrete(name = "") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_pubr() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.position = "right")
