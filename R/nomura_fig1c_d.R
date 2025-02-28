#################################
## Title: Figure 1 panel c-d in Nomura et al - dataset overview
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a cohort-level UMAP and breakdown of the different cell type frequencies
#################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Setup a seurat object for plotting UMAPs
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#############################################################################################################################################################
# Select cell type marker genes
#############################################################################################################################################################

junk_genes <- c(rownames(umi_data_all[[1]])[grep("\\.", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("-AS*", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("LINC", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^RP[S|L]", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^MT-", rownames(umi_data_all[[1]]))])

# Remove uninformative genes (ribosomal proteins, mitochondrially encoded, antisense, LINCs etc.)
selected_genes <- setdiff(rownames(umi_data_all[[1]]), junk_genes)

# Bind all sample-specific matrices to a single one
umi_data <- lapply(unique(meta_data$Sample), function(sname) {
  m <- umi_data_all[[sname]]
  m <- m[selected_genes, meta_data$CellID[meta_data$Sample == sname]]
  m
})
umi_data <- do.call(cbind, umi_data)

# Normalize to counts per million
m <- umi2upm(umi_data)

# Select markers for each cell type
de_genes <- lapply(unique(meta_data$CellType), function(ct) {
  
  print(ct)
  
  d1 <- meta_data %>%
    filter(CellType == ct)
  d2 <- meta_data %>%
    filter(CellType != ct)
  
  m1 <- m[, d1$CellID]
  m2 <- m[, d2$CellID]
  
  rm1 <- log2(rowMeans(m1) + 1)
  rm2 <- log2(rowMeans(m2) + 1)
  
  res <- tibble(CellType = ct, Gene = names(rm1), log2FC = rm1 - rm2)
  
  return(res)
})
de_genes <- do.call(rbind, de_genes)

selected_genes <- de_genes %>%
  filter(Gene %ni% junk_genes) %>%
  filter(CellType %in% c("Malignant", "Astrocyte", "Excitatory neuron", "Inhibitory neuron", "OPC", "Oligodendrocyte", "Endothel", "Pericyte", "Macrophage", "Tcell", "Bcell")) %>%
  group_by(CellType) %>%
  arrange(CellType, desc(log2FC)) %>%
  top_n(100, log2FC)
table(selected_genes$Gene) %>% sort()

selected_genes <- unique(selected_genes$Gene)

#############################################################################################################################################################
# Generate UMI matrix for seurat dataset
#############################################################################################################################################################

# Generate the seurat object with cell type-specific markers
sds <- new_seurat_object(m = umi_data[selected_genes, meta_data$CellID], md = meta_data, project = "GBM-Longitudinal_integrated_analysis_final", verbose = T)

# Re-run umap and clustering
sds <- RunUMAP(sds, dims = 1:100, metric = "correlation", reduction.name = "umap_corr")
sds <- FindNeighbors(sds, dims = 1:2, reduction = "umap_corr", graph.name = "umap_corr_graph")
sds <- FindClusters(sds, resolution = .5, graph.name = "umap_corr_graph")

d <- as_tibble(sds@meta.data)

umap_res <- Reductions(sds, "umap_corr")
umap_res <- Reductions(sds, "umap")
d$UMAP1 <- umap_res@cell.embeddings[d$CellID, 1]
d$UMAP2 <- umap_res@cell.embeddings[d$CellID, 2]

sds <- sds[, meta_data$CellID]

d <- as_tibble(sds@meta.data)

d <- d %>%
  filter(CellID %in% meta_data$CellID)

umap_res <- Reductions(sds, "umap_corr")
d$UMAP1 <- umap_res@cell.embeddings[d$CellID, 1]
d$UMAP2 <- umap_res@cell.embeddings[d$CellID, 2]

# Aggregate unimportant cell types to "Other" category
d <- d %>%
  mutate(CellType = as.character(CellType),
         CellType = case_when(CellType %in% c("Other neuron", "Other normal", "Unresolved") ~ "Other",
                              CellType %in% c("Bcell", "Tcell") ~ "Lymphocyte",
                              CellType == "Macrophage" ~ "TAM",
                              TRUE ~ CellType),
         CellType = factor(CellType, names(celltype_color_vec_reduced)))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 1c - UMAP of all cells colored by cell type
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ggplot(d, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  ggrastr::geom_point_rast(size = .25) +
  scale_color_manual(name = "", values = celltype_color_vec_reduced) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_pubr() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.position = "right")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 1D - Cell type composition pie chart
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

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

d_stats <- d %>% 
  group_by(CellType) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(N = sum(n), Freq = n / N, label = paste0(round(Freq * 100, 1), "%"), label_full = paste0(CellType, "(", label, ")")) %>%
  arrange(Freq)

d_stats$CellType <- factor(d_stats$CellType, d_stats$CellType[order(d_stats$Freq)])

d_stats2 <- d_stats %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

ggplot(d_stats, aes(x = "" , y = Freq, fill = CellType)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(name = "", values = celltype_color_vec_reduced) +
  ggrepel::geom_label_repel(data = d_stats2,
                            aes(y = pos, label = label_full),
                            size = 6, nudge_x = 1, show.legend = FALSE) +
  theme_void() +
  theme(legend.text = element_text(size = 16), legend.position = "none")
