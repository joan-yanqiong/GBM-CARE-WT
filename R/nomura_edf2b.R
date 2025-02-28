#################################
## Title: Extended Data Figure 2 panel b in Nomura et al - % samples contributing per cell type per cluster
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Box plot of % samples contributing per cell type per cluster
#################################

sds <- FindNeighbors(sds, dims = 1:2, reduction = "umap_corr")
sds <- FindClusters(sds, resolution = 5)

d <- as_tibble(sds@meta.data)
table(d$seurat_clusters)

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

n_samples_per_celltype <- d %>%
  group_by(Sample, CellType) %>%
  summarise(n = n()) %>%
  group_by(CellType) %>%
  summarise(N = n(), .groups = "drop")

x <- sapply(levels(d$seurat_clusters), function(x) length(unique(d$Sample[d$seurat_clusters == x])))
x <- tibble(seurat_clusters = names(x), n = x) #%>%

dm <- d %>%
  group_by(seurat_clusters, CellType) %>%
  summarise(n = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  filter(Freq == max(Freq)) %>%
  ungroup() %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  select(seurat_clusters, CellType)

dm <- dm %>%
  left_join(x)

dm <- dm %>%
  left_join(n_samples_per_celltype, by = "CellType")

dm$Freq <- dm$n / dm$N

dm <- dm %>%
  arrange(Freq) %>%
  mutate(seurat_clusters = factor(seurat_clusters, seurat_clusters))

inc_celltypes <- dm %>%
  group_by(CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 5) %>%
  pull(CellType)

median_ord <- dm %>%
  group_by(CellType) %>%
  summarise(med = median(Freq), .groups = "drop") %>%
  arrange(med) %>%
  pull(CellType) %>%
  as.character()

dm %>%
  filter(CellType %in% inc_celltypes, CellType != "Other") %>%
  mutate(CellType = factor(as.character(CellType), median_ord)) %>%
  ggboxplot(x = "CellType", y = "Freq", color = "CellType", add = "jitter") +
  scale_color_manual(name = "", values = celltype_color_vec_reduced) +
  xlab("Cell type") +
  ylab("% samples contributing") +
  scale_y_continuous(labels = percent) +
  # guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_pubr() +
  theme(axis.text.x = element_text(size = 16, angle = 90), axis.title = element_text(size = 20),
        legend.position = "right")
