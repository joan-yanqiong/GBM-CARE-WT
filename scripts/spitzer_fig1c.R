#################################
## Title: Spitzer Figure 1 panel c - UMAP of all cells colored by cell type and timepoint
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a UMAP of all cells colored by cell type and timepoint
#################################

####################################################################################################################################
# Setup data
####################################################################################################################################

library(ggraster)

d <- as_tibble(sds@meta.data)

d <- d %>%
  filter(CellID %in% meta_data$CellID)

umap_res <- Reductions(sds, "umap_corr")
d$UMAP1 <- umap_res@cell.embeddings[d$CellID, 1]
d$UMAP2 <- umap_res@cell.embeddings[d$CellID, 2]

d <- d %>%
  mutate(CellType = as.character(CellType),
         CellType = case_when(CellType %in% c("Other neuron", "Other normal", "Unresolved") ~ "Other",
                              CellType %in% c("Bcell", "Tcell") ~ "Lymphocyte",
                              CellType == "Macrophage" ~ "TAM",
                              TRUE ~ CellType),
         CellType = factor(CellType, names(celltype_color_vec_reduced)))

####################################################################################################################################
# Figure 1c (upper panel) - UMAP of all cells colored by cell type
####################################################################################################################################

ggplot(d, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  ggrastr::geom_point_rast(size = .25) +
  scale_color_manual(name = "", values = celltype_color_vec_reduced) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_pubr() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.position = "right")

####################################################################################################################################
# Figure 1c (lower panel) - UMAP of all cells colored by time point
####################################################################################################################################

ggplot(d, aes(x = UMAP1, y = UMAP2, color = PvsR)) +
  ggrastr::geom_point_rast(size = .25) +
  scale_color_discrete(name = "") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_pubr() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.position = "right")
