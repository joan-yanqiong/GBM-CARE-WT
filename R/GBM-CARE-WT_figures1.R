
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Global definitions
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

celltype_color_vec <- c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Macrophage" = "#80B1D3", # 4, 7, 5
                        "Astrocyte" = "purple", "Excitatory neuron" = "#BC80BD", "Inhibitory neuron" = "#FFED6F", "Other neuron" = "#CCEBC5", # 6, 10, 12
                        "Endothel" = "#FCCDE5", "Pericyte" = "#FFFFB3", "Tcell" = "#8DD3C7", "Bcell" = "#BEBADA", "OPC" = "#FDB462", # 8, 2, 1, 3, 11
                        "Other normal" = "#D9D9D9", "Unresolved" = "grey") # 10

celltype_color_vec_reduced <- c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "TAM" = "#80B1D3",
                                "Astrocyte" = "#BEBADA", "Excitatory neuron" = "#BC80BD", "Inhibitory neuron" = "#FFED6F",
                                "Endothel" = "#FCCDE5", "Pericyte" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "OPC" = "#CCEBC5",
                                "Other" = "grey") # 10

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 1 - Dataset overview
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

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

selected_genes <- setdiff(rownames(umi_data_all[[1]]), junk_genes)

umi_data <- lapply(unique(meta_data$Sample), function(sname) {
  m <- umi_data_all[[sname]]
  m <- m[selected_genes, meta_data$CellID[meta_data$Sample == sname]]
  m
})
umi_data <- do.call(cbind, umi_data)

m <- umi2upm(umi_data)

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

sds <- new_seurat_object(m = umi_data[selected_genes, meta_data$CellID], md = meta_data, project = "GBM-Longitudinal_integrated_analysis_final", verbose = T)

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

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF1a - QC metrics
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- meta_data %>%
  as_tibble()

ggplot(d, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, color = "black", fill = "dodgerblue", alpha = .25) +
  geom_vline(xintercept = quantile(d$nCount_RNA, .1), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = quantile(d$nCount_RNA, .9), linetype = "dashed", color = "black", size = 1) +
  xlab("Total counts per cell [log10]") +
  ylab("Number of cells") +
  scale_x_log10() +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

ggplot(d, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, color = "black", fill = "dodgerblue", alpha = .25) +
  geom_vline(xintercept = quantile(d$nFeature_RNA, .1), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = quantile(d$nFeature_RNA, .9), linetype = "dashed", color = "black", size = 1) +
  xlab("Number of detected genes per cell") +
  ylab("Number of cells") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

ggplot(d, aes(x = Percent.MT)) +
  geom_histogram(bins = 100, color = "black", fill = "dodgerblue", alpha = .25) +
  xlab("Mitochondrial expression [%]") +
  ylab("Number of cells") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF1c - Sample-level CNA matrix
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- cna_sig_per_chr %>%
  filter(CellID %in% mdata$CellID)

d <- d %>%
  group_by(Patient, Timepoint, Sample, CHR) %>%
  summarise(Sig = mean(Sig))
d$ID <- paste0(d$Patient, d$Timepoint)

d$Sig[abs(d$Sig) < .05] <- 0

m <- acast(d, ID ~ CHR, value.var = "Sig")

hc <- fastcluster::hclust(d = as.dist(1 - cor(t(m))), method = "average")

d$ID <- factor(as.character(d$ID), hc$labels[hc$order])

ggplot(d, aes(x = CHR, y = ID, fill = Sig)) +
  facet_grid(cols = vars(CHR), scales = "free_x") +
  geom_tile() +
  scale_fill_gradient2(name = "", low = "dodgerblue", mid = "white", high = "red",
                       labels = c("DEL", "", "", "", "AMP"), limits = c(-.2, .2), oob = squish) +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1, ticks.colour = "black")) +
  xlab("Chromosomes") +
  ylab("Samples") +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 legend.text = element_text(size = 16),
                 strip.text = element_text(size = 10))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF2a - UMAP of all cells colored by data generating lab
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ggplot(d, aes(x = UMAP1, y = UMAP2, color = Lab)) +
  ggrastr::geom_point_rast(size = .25) +
  scale_color_discrete(name = "") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_pubr() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.position = "right")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF2b - Box plot of % samples contributing per cell type per cluster
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

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

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF2c - Heatmap of cell type markers
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ct_markers <- de_genes %>%
  filter(CellType %ni% c("Other neuron", "Other normal", "Unresolved")) %>%
  arrange(desc(log2FC)) %>%
  group_by(CellType) %>%
  top_n(50, log2FC) %>%
  ungroup()

md <- meta_data %>%
  filter(CellType %in% ct_markers$CellType) %>%
  as_tibble()

umi_data <- lapply(unique(md$Sample), function(sname) {
  print(sname)
  d <- md %>%
    filter(Sample == sname)
  x <- umi_data_all[[sname]][, d$CellID]
  x <- umi2upm(x)
  x <- x[unique(ct_markers$Gene), ]
  x <- log2(x / 10 + 1)
  x
})
umi_data <- do.call(cbind, umi_data)

umi_data <- umi_data[, md$CellID]

tmp <- lapply(unique(ct_markers$CellType), function(ct) {
  print(ct)
  sig <- ct_markers$Gene[ct_markers$CellType == ct]
  colMeans(umi_data[sig, ])
})
tmp <- do.call(rbind, tmp)
rownames(tmp) <- unique(ct_markers$CellType)
tmp <- tmp[, md$CellID]

md$Score <- sapply(1:nrow(md), function(i) {
  cid <- md$CellID[i]
  ct <- as.character(md$CellType[i])
  tmp[ct, cid]
})

d <- md
d <- lapply(unique(d$CellType), function(ct) {
  print(ct)
  x <- d %>%
    filter(CellType == ct)
  if(nrow(x) < 100) {
    x <- x %>%
      arrange(desc(Score))
  } else if(nrow(x) <= 1000) {
    x <- x %>%
      arrange(desc(Score)) %>%
      head(100)
  } else {
    x <- x %>%
      arrange(desc(Score)) %>%
      head(1000)
  }
  x
})
d <- do.call(rbind, d)

m <- umi_data[ct_markers$Gene %>% unique(),
              d$CellID]

m <- rowcenter(m)

dm <- melt(m) %>%
  as_tibble() %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
colnames(dm) <- c("Gene", "CellID", "value")

dm  <- dm %>%
  left_join(d %>%
              select(CellID, CellType))
colnames(dm)[ncol(dm)] <- "CellType.x"

dm  <- dm %>%
  left_join(ct_markers %>%
              select(Gene, CellType) %>%
              filter(!duplicated(Gene)))
colnames(dm)[ncol(dm)] <- "CellType.y"

marker_genes <- rbind(tibble(CellType = "Malignant", Gene = c("EGFR", "MEOX2", "PTPRZ1", "ETV1", "METTL7B")),
                      tibble(CellType = "Oligodendrocyte", Gene = c("MOG", "MOBP", "MAG", "OPALIN", "PLP1")),
                      tibble(CellType = "Macrophage", Gene = c("CD163", "MSR1", "CSF2RA", "FCGR2A", "IL1B")),
                      tibble(CellType = "Astrocyte", Gene = c("SLC1A2", "NRG3", "MAOB", "FGFR3", "ADGRV1")),
                      tibble(CellType = "Excitatory neuron", Gene = c("SYT1", "RBFOX3", "GRIN1", "GRIN2A", "GRM1")),
                      tibble(CellType = "Inhibitory neuron", Gene = c("GAD1", "GAD2", "VIP", "SST", "GRIK1")),
                      tibble(CellType = "Endothel", Gene = c("FLT1", "VWF", "CLDN5", "ESM1", "ESAM")),
                      tibble(CellType = "Pericyte", Gene = c("ENPEP", "COL1A1", "PDGFRB", "ACTA2", "RGS5")),
                      tibble(CellType = "Tcell", Gene = c("CD2", "GZMA", "IL7R", "CD3E", "ZAP70")),
                      tibble(CellType = "Bcell", Gene = c("IGHM", "IGKC", "IGHG1", "MS4A1", "JCHAIN")),
                      tibble(CellType = "OPC", Gene = c("PDGFRA", "FGF14", "OLIG2", "SOX11", "MYT1")))
marker_genes$CellType <- factor(marker_genes$CellType, levels(de_genes$CellType))

dm <- de_genes %>%
  filter(Gene %in% marker_genes$Gene, CellType %ni% c("Other neuron", "Other normal", "Unresolved"))

dm <- dm %>%
  left_join(marker_genes, by = "Gene")

ggplot(data = dm,
       aes(x = CellType.x, y = Gene, fill = log2FC)) + 
  geom_tile(color = "black") + 
  facet_grid(rows = vars(CellType.y), scales = "free", space = "free_x") +
  scale_fill_gradient2(name = "Relative expression [log2]", low = "dodgerblue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-4, 4), oob = squish) +
  xlab("Cell types") +
  ylab("Marker genes") +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1, ticks.colour = "black")) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 12, angle = 90),
                 axis.text.y = element_text(size = 10), strip.text = element_text(size = 12))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2 - Revisiting malignant states
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2a - MP heatmap
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

dm <- malignant_nmf_metaprograms$nmf_intersect_meltI %>%
  as_tibble()
clusters <- lapply(names(malignant_nmf_metaprograms$clusters), function(xn) {
  tibble(Var1 = malignant_nmf_metaprograms$clusters[[xn]], Cluster = xn)
})
clusters <- do.call(rbind, clusters)
dm$Var1 <- as.character(dm$Var1)
dm <- dm %>%
  left_join(clusters, by = "Var1")
dm$Var1 <- factor(dm$Var1, levels(dm$Var2))

exc_p <- dm %>%
  filter(Cluster %in% c("Cluster_1", "Cluster_11", "Cluster_12")) %>%
  pull(Var1) %>%
  as.character() %>%
  unique()

dm <- dm %>%
  filter(as.character(Var1) %ni% exc_p, as.character(Var2) %ni% exc_p)

dm$Var1 <- droplevels(dm$Var1)
dm$Var2 <- droplevels(dm$Var2)

limits <- lengths(malignant_nmf_metaprograms$clusters)
limits <- limits[names(limits) %ni% c("Cluster_1", "Cluster_11", "Cluster_12")]
limits <- limits %>% cumsum()
limits <- c(1, limits[-length(limits)])

x1 <- limits[seq(from = 1, to = length(limits) - 1, by = 1)]
x2 <- limits[seq(from = 2, to = length(limits), by = 1)]
x1 <- factor(levels(dm$Var1)[x1], levels(dm$Var1))
x2 <- factor(levels(dm$Var1)[x2], levels(dm$Var1))
y1 <- x1; y2 <- x2
f <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)

ggplot(data = dm,
       aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  geom_rect(data = f, aes(x = NULL,y = NULL, xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            color = "black", linetype = "dashed", fill = NA, size = .75, inherit.aes = F) +
  scale_y_discrete(limits = rev)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2b - MP scores heatmap
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- mdata

sampled_cells <- d %>% 
  sample_n(20000)

sampled_cells <- sampled_cells %>%
  dplyr::select(CellID, Sample, Patient, Timepoint, starts_with("MP"), State)

sampled_cells <- sampled_cells %>%
  rename(Cilia = MP_13_Cilia,
         AC = MP_4_AC,
         Hypoxia = MP_5_Hypoxia,
         MES = MP_6_MES,
         GPC = MP_8_GPC,
         OPC = MP_2_OPC,
         NPC = MP_7_NPC,
         CC = MP_3_CC)

sampled_cells$Neuron <- sapply(1:nrow(sampled_cells), function(i) max(sampled_cells[i, c("MP_9_ExN", "MP_14_NRGN")]))
sampled_cells$Stress <- sapply(1:nrow(sampled_cells), function(i) max(sampled_cells[i, c("MP_10_Stress1", "MP_15_Stress2")]))

sampled_cells <- sampled_cells %>%
  select(-starts_with("MP_"))

scores <- sampled_cells %>% dplyr::select(Cilia, AC, MES, Hypoxia, Stress, GPC, OPC, NPC, Neuron, CC)
scores <- as.matrix(scores)
rownames(scores) <- sampled_cells$CellID
scores <- t(scores)

dm <- melt(scores) %>%
  as_tibble()
colnames(dm) <- c("MP", "CellID", "Score")

dm <- dm %>%
  left_join(sampled_cells %>%
              select(CellID, State), by = "CellID")

dm$State <- factor(dm$State, c("Cilia", "AC", "MES", "Hypoxia", "Stress", "GPC", "OPC", "NPC", "Neuron", "Unresolved"))

dm$MP <- factor(dm$MP, c("Cilia", "AC", "MES", "Hypoxia", "Stress", "GPC", "OPC", "NPC", "Neuron", "CC"))

dm$CellID <- factor(dm$CellID, colnames(scores)[order(scores["CC", ])])

p1 <- dm %>%
  filter(MP != "CC") %>%
  ggplot(aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(State), scales = "free", space = "free_x") +
  geom_tile() +
  scale_fill_gradient2(name = "Score", low = "dodgerblue", mid = "white", high = "red",
                       labels = c("-1", "", "0", "", "1"), limits = c(-1, 1), oob = squish) +
  xlab("") +
  ylab("States") +
  guides(fill = guide_colourbar(ticks.colour = "black", frame.colour = "black")) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 12))

p2 <- dm %>%
  filter(MP == "CC") %>%
  ggplot(aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(State), scales = "free", space = "free_x") +
  geom_tile() +
  geom_raster() +
  scale_fill_gradient2(name = "Score", low = "dodgerblue", mid = "white", high = "red",
                       labels = c("-1", "", "0", "", "1"), limits = c(-1, 1), oob = squish) +
  xlab("Cells") +
  ylab("States") +
  guides(fill = guide_colourbar(ticks.colour = "black", frame.colour = "black")) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 strip.text = element_text(size = 12), legend.position = "none")

p1 + p2 + plot_layout(nrow = 2, heights = c(5, .5))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2c - Heatmap of new states
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

new_states <- mdata %>%
  filter(State_full %in% c("MP_8_GPC", "MP_9_ExN", "MP_13_Cilia"))

new_states_scores <- state_data %>%
  ungroup() %>%
  filter(CellID %in% new_states$CellID, Program %in% c("MP_8_GPC", "MP_9_ExN", "MP_13_Cilia"), p.sig == T)

is_dup <- new_states_scores %>%
  group_by(CellID) %>%
  filter(duplicated(CellID)) %>%
  pull(CellID)

new_states_scores <- new_states_scores %>%
  filter(CellID %ni% is_dup)

new_states_scores$State <- new_states_scores$Class

selected_cells <- new_states_scores %>%
  group_by(Program) %>%
  arrange(desc(Score)) %>%
  top_frac(.05, Score) %>%
  ungroup()

reference_cells <- mdata %>%
  filter(State %in% c("AC", "MES", "Hypoxia", "OPC", "NPC")) %>%
  group_by(State) %>%
  sample_n(1000)

cells <- rbind(selected_cells %>%
                 select(CellID, Sample),
               reference_cells %>%
                 ungroup() %>%
                 select(CellID, Sample))

mp_genes <- rbind(tibble(Gene = MP_list_named$MP_13_Cilia, MP = "Cilia"),
                  tibble(Gene = MP_list_named$MP_8_GPC, MP = "GPC"),
                  tibble(Gene = MP_list_named$MP_9_ExN, MP = "Neuron"))

m <- lapply(unique(cells$Sample), function(sname) {
  print(sname)
  res <- umi_data_all[[sname]][mp_genes$Gene[mp_genes$Gene %in% rownames(umi_data_all[[sname]])], cells$CellID[cells$Sample == sname]]
  if(is.null(dim(res))) {
    res <- as.matrix(res)
    colnames(res) <- cells$CellID[cells$Sample == sname]
    rownames(res) <- mp_genes$Gene[mp_genes$Gene %in% rownames(umi_data_all[[sname]])]
  }
  return(res)
})
m <- do.call(cbind, m)

m <- umi2upm(m)
m <- log2(m / 10 + 1)
m <- rowcenter(m)

m <- m[, selected_cells$CellID]

dm <- melt(m) %>%
  as_tibble()
colnames(dm) <- c("Gene", "CellID", "value")

dm <- dm %>%
  left_join(selected_cells %>%
              select(CellID, State),
            by = "CellID") %>%
  mutate(State = factor(State, c("Cilia", "GPC", "Neuron")))

dm <- dm %>%
  left_join(mp_genes %>%
              select(Gene, MP),
            by = "Gene") %>%
  mutate(MP = factor(MP, c("Cilia", "GPC", "Neuron")))

ggplot(dm, aes(x = CellID, y = Gene, fill = value)) +
  facet_grid(rows = vars(MP), cols = vars(State), scales = "free") +
  geom_tile() +
  scale_fill_distiller(name = "Centered log2-expression", palette = "RdBu", direction = -1,
                       breaks = seq(-5, 5),
                       limits = c(-5, 5), oob = squish, guide = "colorsteps") +
  xlab("Cells") +
  ylab("Genes") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2d - State similarity with states from normal brain development
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

LIU_markers <- read.csv("hNSPC_marker_genes.csv", header = T, stringsAsFactors = F)
LIU_markers <- as_tibble(LIU_markers)

LIU_markers <- apply(LIU_markers, 2, function(x) {
  sapply(strsplit(x, split = "_"), function(y) y[[2]][1]) %>%
    head(50)
})

LIU_markers_list <- setNames(lapply(1:ncol(LIU_markers), function(i) LIU_markers[, i]), colnames(LIU_markers))

LIU_scores <- score_within_samples(umi_data_list = umi_data_all, md = mdata, sigs = LIU_markers_list)

m1 <- MP_scores[, -c(1:4)]
m1 <- m1 %>%
  select(starts_with("MP_"), -MP_1_RP, -MP_11_MIC, -MP_12_LQ, -MP_16_GlioNeural)
m1 <- as.matrix(m1)
rownames(m1) <- MP_scores$CellID

m2 <- LIU_scores[, -c(1:4)]
m2 <- as.matrix(m2)
rownames(m2) <- LIU_scores$CellID

m1 <- m1[rownames(m2), ]

cor_m1 <- cor(m1, m2)

load(file = "BrainSignatures.rda")

NBD_scores <- score_within_samples(umi_data_list = umi_data_all, md = mdata, sigs = BrainSignatures)

m1 <- MP_scores[, -c(1:4)]
m1 <- m1 %>%
  select(starts_with("MP_"), -MP_1_RP, -MP_11_MIC, -MP_12_LQ, -MP_16_GlioNeural)
m1 <- as.matrix(m1)
rownames(m1) <- MP_scores$CellID

m2 <- NBD_scores[, -c(1:4)]
m2 <- as.matrix(m2)
rownames(m2) <- NBD_scores$CellID

m1 <- m1[rownames(m2), ]

cor_m2 <- cor(m1, m2)

mps <- c("MP_6_MES", "MP_13_Cilia", "MP_4_AC", "MP_8_GPC", "MP_2_OPC", "MP_7_NPC", "MP_9_ExN", "MP_3_CC")

liu <- c("Ventricular.radial.glia", "Outer.radial.glia",
         "Glial.progenitor.cell", "Astrocyte",
         "OPC..dividing.", "Pre.OPC", "OPC", "Oligodendrocyte",
         "Intermediate.progenitor", "Excitatory.neuron..early.", "Excitatory.neuron..late.", "Inhibitory.neuron", "Subcortical.neuron",
         "Radial.glia..S.phase.", "Radial.glia..G2M.phase.")

nbd <- c("NEU.EX.L2.3-velm19", "Choroid-nowak17", "RG.v-nowak17", "RG.o-nowak17",
         "RG.o-poli19", "RG.v-poli19", "AC.FB-velm19")

liu_nbd_ord <- c("Choroid-nowak17", "Astrocyte", "AC.FB-velm19",
                 "RG.v-nowak17", "RG.v-poli19", "Ventricular.radial.glia",
                 "RG.o-poli19", "RG.o-nowak17", "Outer.radial.glia",
                 "Glial.progenitor.cell",
                 "OPC..dividing.", "Pre.OPC", "OPC", "Oligodendrocyte",
                 "Intermediate.progenitor", "Excitatory.neuron..early.", "Excitatory.neuron..late.", "Inhibitory.neuron", "Subcortical.neuron",
                 "NEU.EX.L2.3-velm19",
                 "Radial.glia..S.phase.", "Radial.glia..G2M.phase.")

cor_m <- cbind(cor_m1[mps, liu],
               cor_m2[mps, nbd])
dm <- melt(cor_m)

dm <- dm[dm$Var1 %in% mps, ]
dm$Var1 <- factor(dm$Var1, mps)

dm <- dm[dm$Var2 %in% liu_nbd_ord, ]
dm$Var2 <- factor(dm$Var2, liu_nbd_ord)

ggplot(dm, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Correlation", low = "dodgerblue", mid = "white", high = "red", limits = c(-.75, .75), oob = squish,
                       breaks = c(-.75, -.5, -.25, 0, .25, .5, .75), labels = c("-.75", "", "", "0", "", "", ".75")) +
  xlab("Normal brain states") +
  ylab("GBM states") +
  scale_x_discrete(labels = c("Choroid (Nowakowski)",
                              "Astrocyte (Liu)", "Astrocyte (Velmeshev)",
                              "vRG (Nowakowski)", "vRG (Polioudakis)", "vRG (Liu)",
                              "oRG (Nowakowski)", "oRG (Polioudakis)", "oRG (Liu)",
                              "GPC (Liu)",
                              "Dividing OPC (Liu)", "Pre-OPC (Liu)", "OPC (Liu)", "Oligodendrocyte (Liu)",
                              "IPC (Liu)", "Ex. neuron - early (Liu).",
                              "Ex. neuron - late (Liu)", "In. neuron (Liu)", "Subcortical neuron (Liu)",
                              "Ex. neuron (Velmeshev)",
                              "Cell cycle - G1S (Liu)", "Cell cycle - G2M (Liu)")) +
  scale_y_discrete(labels = c("MES-like", "Cilia-like", "AC-like", "GPC-like", "OPC-like",
                              "NPC-like", "NEU-like", "Cell cycle")) +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 12, angle = 90))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2f - Example of hybrid states (AC, GPC, OPC)
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
# Select the singular cells
####################################################################################################################################

singular_states <- state_data %>%
  ungroup() %>%
  filter(Program %in% c("MP_4_AC", "MP_8_GPC"), p.sig == T) %>%
  filter(CellID %ni% state_data_hybrid$CellID)

singular_states <- singular_states %>%
  group_by(Program) %>%
  arrange(Program, desc(Score)) %>%
  top_n(100, Score) %>%
  ungroup()

####################################################################################################################################
# Select the AC-GPC hybrids
####################################################################################################################################

ac_gpc_hybrids <- state_data %>%
  ungroup() %>%
  filter(Program %in% c("MP_4_AC", "MP_8_GPC"), p.sig == T) %>%
  group_by(CellID) %>%
  filter("MP_4_AC" %in% Program & "MP_8_GPC" %in% Program) %>%
  ungroup() %>%
  filter(CellID %in% state_data_hybrid$CellID[state_data_hybrid$N == 2])

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  filter(Program == "MP_4_AC") %>%
  rename(AC = Score) %>%
  select(CellID, AC) %>%
  left_join(ac_gpc_hybrids %>%
              filter(Program == "MP_8_GPC") %>%
              rename(GPC = Score) %>%
              select(CellID, GPC),
            by = "CellID")
ac_gpc_hybrids$Diff <- abs(ac_gpc_hybrids$AC - ac_gpc_hybrids$GPC)

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  mutate(Sum = AC + GPC)

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  arrange(Diff)

gghistogram(data = ac_gpc_hybrids, x = "Diff", bins = 50)

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  filter(Diff < .01) %>%
  arrange(desc(Sum)) %>%
  top_n(100, Diff)

####################################################################################################################################
# Generate the heatmap
####################################################################################################################################

hybrid_vs_singular_cells <- rbind(singular_states %>%
                                    select(CellID, Program) %>%
                                    mutate(State = case_when(Program == "MP_4_AC" ~ "AC",
                                                             Program == "MP_8_GPC" ~ "GPC")) %>%
                                    select(CellID, State),
                                  ac_gpc_hybrids %>%
                                    select(CellID) %>%
                                    mutate(State = "AC-GPC"))

hybrid_vs_singular_cells$Sample <- get_sname(hybrid_vs_singular_cells$CellID, split = "-")

mp_genes <- rbind(tibble(Gene = MP_list_named$MP_8_GPC, MP = "GPC"),
                  tibble(Gene = MP_list_named$MP_4_AC, MP = "AC"))

m <- lapply(unique(hybrid_vs_singular_cells$Sample), function(sname) {
  print(sname)
  res <- umi_data_all[[sname]][mp_genes$Gene, hybrid_vs_singular_cells$CellID[hybrid_vs_singular_cells$Sample == sname]]
  if(is.null(dim(res))) {
    res <- as.matrix(res)
    colnames(res) <- hybrid_vs_singular_cells$CellID[hybrid_vs_singular_cells$Sample == sname]
    rownames(res) <- mp_genes$Gene
  }
  return(res)
})
m <- do.call(cbind, m)

m <- umi2upm(m)
m <- log2(m / 10 + 1)
m <- rowcenter(m)

dm <- melt(m) %>%
  as_tibble()
colnames(dm) <- c("Gene", "CellID", "value")

dm <- dm %>%
  left_join(hybrid_vs_singular_cells %>%
              select(CellID, State),
            by = "CellID") %>%
  mutate(State = factor(State, c("AC", "AC-GPC", "GPC")))

dm <- dm %>%
  left_join(mp_genes %>%
              select(Gene, MP),
            by = "Gene") %>%
  mutate(MP = factor(MP, c("AC", "GPC")))

ord1 <- ac_gpc_hybrids %>%
  mutate(Diff = GPC - AC) %>%
  arrange(Diff) %>%
  pull(CellID)

ord2 <- singular_states %>%
  pull(CellID)

dm$CellID <- factor(dm$CellID, c(ord1, ord2))

ggplot(dm, aes(x = CellID, y = Gene, fill = value)) +
  facet_grid(rows = vars(MP), cols = vars(State), scales = "free") +
  geom_tile() +
  scale_fill_distiller(name = "Centered log2-expression", palette = "RdBu", direction = -1,
                       breaks = seq(-3, 3),
                       limits = c(-3, 3), oob = squish, guide = "colorsteps") +
  xlab("Cells") +
  ylab("Genes") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2g - Lineage model inferred from hybrdis
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- obs_vs_exp_hybrids %>%
  filter(sigRes)

edges <- lapply(d$Hybrid, function(h) {
  
  h <- strsplit(x = as.character(h), split = "-")
  
  vs1 <- h[[1]][1]
  vs2 <- h[[1]][2]
  data.frame(S1 = vs1, S2 = vs2)
})
edges <- do.call(rbind, edges)
edges$Freq <- d$ObsFreq
edges$log2FC <- d$log2FC_2
edges$FC <- 2^edges$log2FC

edges <- edges %>%
  filter(S2 != "Stress")

hybrid_graph <- graph_from_data_frame(edges, directed = F)

E(hybrid_graph)$weight <- edges$log2FC
E(hybrid_graph)$arrow.size <- edges$log2FC

plot(hybrid_graph,
     edge.arrow.size = E(hybrid_graph)$arrow.size, edge.width = E(hybrid_graph)$weight, edge.label.cex = 2,
     edge.label = E(hybrid_graph)$label, edge.curved = F,
     vertex.size = 25, vertex.label.cex = 2)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 2h - GPC butterfly
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

valid_patients <- mdata %>%
  group_by(Patient, Timepoint) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n >= 100) %>%
  filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
  pull(Patient) %>%
  unique()

d <- mdata %>%
  filter(Patient %in% valid_patients) %>%
  group_by(Patient, Timepoint) %>%
  sample_n(100)

knn_res <- FNN::get.knn(d[, c("Dx", "Dy")], k = 10)
d$isGPC <- d$State == "GPC"
d$GPC_density <- sapply(1:nrow(knn_res$nn.index), function(i) sum(d$isGPC[knn_res$nn.index[i, ]]))
table(d$GPC_density)

d$GPC_density <- d$GPC_density / 10

d <- d %>% arrange(GPC_density)

ggplot(data = d, mapping = aes(x = Dx, y = Dy)) +
  geom_point(mapping = aes(fill = GPC_density), shape = 21, color = "black", size = 4,
             show.legend = c("color" = TRUE, "fill" = TRUE, "alpha" = FALSE)) +
  scale_fill_gradient2(name = "% GPC neighbors", low = "black", mid = "brown", high = "red", oob = squish, limits = c(0, .4), midpoint = .2, labels = c("0%", "", "20%", "", "40%")) +
  xlab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
  ylab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr() +
  theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF3a - MP heatmap with low-quality MPs
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

limits <- lengths(malignant_nmf_metaprograms$clusters) %>% cumsum()
limits <- c(1, limits[-length(limits)])

x1 <- limits[seq(from = 1, to = length(limits) - 1, by = 1)]
x2 <- limits[seq(from = 2, to = length(limits), by = 1)]
x1 <- factor(levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1)[x1],
             levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1))
x2 <- factor(levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1)[x2],
             levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1))
y1 <- x1; y2 <- x2
f <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)

ggplot(data = malignant_nmf_metaprograms$nmf_intersect_meltI,
       aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  geom_rect(data = f, aes(x = NULL,y = NULL, xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            color = "black", linetype = "dashed", fill = NA, size = .75, inherit.aes = F) +
  scale_y_discrete(limits = rev)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF3b-d - Comparison with Neftel data
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
# Prepare data for plotting
####################################################################################################################################

nl_exp_data <- read.delim("IDHwtGBM.processed.SS2.logTPM.txt", header = T)

rownames(nl_exp_data) <- nl_exp_data[, 1]
nl_exp_data <- nl_exp_data[, -1]
nl_exp_data <- as.matrix(nl_exp_data)
colnames(nl_exp_data) <- gsub("\\.", "-", colnames(nl_exp_data))

nl_meta_data <- read.delim("IDHwt.GBM.Metadata.SS2.txt", header = T, stringsAsFactors = F)

nl_meta_data <- nl_meta_data[-1, ]
nl_meta_data <- as_tibble(nl_meta_data)

nl_meta_data <- nl_meta_data %>%
  filter(CellAssignment == "Malignant")

sigs <- c(setNames(Signatures_GBM, paste0("NL_", names(Signatures_GBM))),
          MP_list_named)

nl_scores <- lapply(unique(nl_meta_data$Sample), function(s_name) {
  
  d <- nl_meta_data %>% filter(Sample == s_name)
  
  x <- nl_exp_data[, d$NAME]
  dim(x)
  
  res <- scalop::sigScores(m = x, sigs = sigs, conserved.genes = .25)
  res <- as_tibble(res, rownames = "NAME")
  
  res$Sample <- s_name
  
  return(res)
})
nl_scores <- do.call(rbind, nl_scores)

nl_states <- c("NL_AC", "NL_MES1", "NL_MES2", "NL_OPC", "NL_NPC1", "NL_NPC2")
mp_states <- c("MP_2_OPC", "MP_4_AC", "MP_5_Hypoxia", "MP_6_MES", "MP_7_NPC", "MP_8_GPC", "MP_9_ExN", "MP_13_Cilia")

nl_scores$MaxNLscore <- apply(nl_scores[, nl_states], 1, max)
nl_scores$MaxNLstate <- apply(nl_scores[, nl_states], 1, function(x) nl_states[which.max(x)])

nl_scores$MaxMPscore <- apply(nl_scores[, mp_states], 1, max)
nl_scores$MaxMPstate <- apply(nl_scores[, mp_states], 1, function(x) mp_states[which.max(x)])

nl_scores_care <- score_within_samples(umi_data_list = umi_data_all, md = mdata %>% filter(Sample %in% pt_pairs$Sample), sigs = setNames(Signatures_GBM, paste0("NL_", names(Signatures_GBM))))

nl_scores_care <- nl_scores_care %>%
  left_join(MP_scores %>%
              select(CellID, starts_with("MP_")), by = "CellID")

nl_scores_care$MaxNLscore <- apply(nl_scores_care[, nl_states], 1, max)
nl_scores_care$MaxNLstate <- apply(nl_scores_care[, nl_states], 1, function(x) nl_states[which.max(x)])

nl_scores_care$MaxMPscore <- apply(nl_scores_care[, mp_states], 1, max)
nl_scores_care$MaxMPstate <- apply(nl_scores_care[, mp_states], 1, function(x) mp_states[which.max(x)])

####################################################################################################################################
# Figure EDF3b - NL vs. CARE scores cor
####################################################################################################################################

cor_m <- cor(nl_scores_care %>%
               select(starts_with("NL_")),
             nl_scores_care %>%
               select(all_of(mps)))

hc1 <- fastcluster::hclust(d = dist(cor_m, "euclidean"))
hc2 <- fastcluster::hclust(d = dist(t(cor_m), "euclidean"))

cor_m <- melt(cor_m) %>%
  as_tibble()

cor_m$Var1 <- factor(as.character(cor_m$Var1), hc1$labels[hc1$order])
cor_m$Var2 <- factor(as.character(cor_m$Var2), hc2$labels[hc2$order])

ggplot(cor_m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Correlation", low = "dodgerblue", mid = "white", high = "red", limits = c(-.75, .75), oob = squish,
                       breaks = c(-.75, -.5, -.25, 0, .25, .5, .75), labels = c("-.75", "", "", "0", "", "", ".75")) +
  xlab("CARE MPs") +
  ylab("Neftel MPs") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 20, angle = 90))

####################################################################################################################################
# Figure EDF3c - NL vs. CARE Jaccard
####################################################################################################################################

mps <- c("MP_2_OPC", "MP_3_CC", "MP_4_AC", "MP_5_Hypoxia", "MP_6_MES", "MP_7_NPC", "MP_8_GPC", "MP_9_ExN", "MP_10_Stress1", "MP_13_Cilia", "MP_14_NRGN", "MP_15_Stress2")

jac_m <- jaccard(MP_list_named[mps], Signatures_GBM)

hc1 <- fastcluster::hclust(d = dist(jac_m, "euclidean"))
hc2 <- fastcluster::hclust(d = dist(t(jac_m), "euclidean"))

jac_m <- melt(jac_m) %>%
  as_tibble()

jac_m$Var1 <- factor(as.character(jac_m$Var1), hc1$labels[hc1$order])
jac_m$Var2 <- factor(as.character(jac_m$Var2), hc2$labels[hc2$order])

ggplot(jac_m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(name = "Jaccard", low = "white", high = "red", limits = c(0, .15), oob = squish) +
  xlab("CARE MPs") +
  ylab("Neftel MPs") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 20, angle = 90))

####################################################################################################################################
# Figure EDF3d - Compare NL scores with CARE scores
####################################################################################################################################

res1 <- rbind(nl_scores %>%
                filter(MaxNLstate == "NL_AC" & MaxMPstate == "MP_4_AC") %>%
                mutate(State = "AC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_MES1" & MaxMPstate == "MP_6_MES") %>%
                mutate(State = "MES") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_MES2" & MaxMPstate == "MP_5_Hypoxia") %>%
                mutate(State = "Hypoxia") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_OPC" & MaxMPstate == "MP_2_OPC") %>%
                mutate(State = "OPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_NPC2" & MaxMPstate == "MP_7_NPC") %>%
                mutate(State = "NPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt()) %>%
  as_tibble() %>%
  mutate(Dataset = "Neftel")

res2 <- rbind(nl_scores_care %>%
                filter(MaxNLstate == "NL_AC" & MaxMPstate == "MP_4_AC") %>%
                mutate(State = "AC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_MES1" & MaxMPstate == "MP_6_MES") %>%
                mutate(State = "MES") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_MES2" & MaxMPstate == "MP_5_Hypoxia") %>%
                mutate(State = "Hypoxia") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_OPC" & MaxMPstate == "MP_2_OPC") %>%
                mutate(State = "OPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_NPC2" & MaxMPstate == "MP_7_NPC") %>%
                mutate(State = "NPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt()) %>%
  as_tibble() %>%
  mutate(Dataset = "CARE")

res <- rbind(res1, res2)

res_stats <- res %>%
  group_by(Dataset, State, variable) %>%
  summarise(Mean = mean(value), .groups = "drop")

p1 <- res %>%
  filter(Dataset == "CARE") %>%
  ggplot(aes(x = value, y = after_stat(ndensity), fill = variable)) +
  geom_density(color = "black", alpha = .3) +
  geom_vline(data = res_stats %>%
               filter(Dataset == "CARE"),
             mapping = aes(xintercept = Mean, color = variable), linetype = "dashed", size = 1, show.legend = F) +
  facet_grid(rows = vars(Dataset), cols = vars(State)) +
  scale_fill_discrete(name = "Signature origin") +
  xlab("Score") +
  ylab("Density (scaled to 1)") +
  theme_gbm_pvsr()

p2 <- res %>%
  filter(Dataset == "Neftel") %>%
  ggplot(aes(x = value, y = after_stat(ndensity), fill = variable)) +
  geom_density(color = "black", alpha = .3) +
  geom_vline(data = res_stats %>%
               filter(Dataset == "Neftel"),
             mapping = aes(xintercept = Mean, color = variable), linetype = "dashed", size = 1, show.legend = F) +
  facet_grid(rows = vars(Dataset), cols = vars(State)) +
  scale_fill_discrete(name = "Signature origin") +
  xlab("Score") +
  ylab("Density (scaled to 1)") +
  theme_gbm_pvsr()

p1 + p2 + plot_layout(nrow = 2)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF3e - MP heatmap with low-quality MPs
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
# Bhaduri 2020
####################################################################################################################################

krieg_sigs_tbl <- read.csv(paste0(TABLES_ROOT, "kriegstein_2020_supp_table_2.csv")) %>%
  as_tibble()
krieg_cluster_int_tbl <- read.csv(paste0(TABLES_ROOT, "kriegstein_2020_supp_table_2_cluster_int.csv")) %>%
  as_tibble()
krieg_cluster_int_tbl$CellType2 <- paste0(krieg_cluster_int_tbl$CellType, "_", krieg_cluster_int_tbl$cluster)

krieg_sigs_tbl <- krieg_sigs_tbl %>%
  left_join(krieg_cluster_int_tbl)

krieg_sigs_tbl$CellType2 <- paste0(krieg_sigs_tbl$CellType, "_", krieg_sigs_tbl$cluster)

krieg_sigs <- lapply(unique(krieg_cluster_int_tbl$CellType2), function(ct) {
  x <- krieg_sigs_tbl %>%
    filter(CellType2 == ct)
  n_clust <- length(unique(x$cluster))  
  
  x <- x %>%
    group_by(gene) %>%
    summarise(avg_logFC = mean(avg_logFC), n = n(), freq = n / n_clust) %>%
    arrange(desc(n), desc(avg_logFC))
  
  x <- x %>%
    head(50) %>%
    pull(gene)
})
names(krieg_sigs) <- unique(krieg_cluster_int_tbl$CellType2)
krieg_sigs <- krieg_sigs[sort(names(krieg_sigs))]

krieg_sigs <- krieg_sigs[names(krieg_sigs) %ni% c("B Cells", "Dividing B Cells", "Endothelial", "Microglia", "Pericyte",
                                                  "Red blood cells", "Tumor Associated Macrophage", "Unknown")]

KRIEG_scores <- score_within_samples(umi_data_list = umi_data_all, md = mdata, sigs = krieg_sigs)

#############################################################################################################################################################
# Mathur 2024
#############################################################################################################################################################

costello_sigs_tbl <- read.csv(paste0(TABLES_ROOT, "costello_2024_cell_type_sigs.csv")) %>%
  as_tibble()
costello_cluster_int_tbl <- read.csv(paste0(TABLES_ROOT, "cotello_2024_dict.csv")) %>%
  as_tibble()

costello_sigs_tbl <- costello_sigs_tbl %>%
  left_join(costello_cluster_int_tbl %>%
              select(Module, CellType))

costello_sigs <- lapply(unique(costello_cluster_int_tbl$CellType), function(ct) {
  x <- costello_sigs_tbl %>%
    filter(CellType == ct)
  
  x <- x %>%
    # head(50) %>%
    pull(Gene)
})
names(costello_sigs) <- unique(costello_cluster_int_tbl$CellType)
costello_sigs <- costello_sigs[sort(names(costello_sigs))]

COSTELLO_scores <- score_within_samples(umi_data_list = umi_data_all, md = mdata, sigs = costello_sigs)

#############################################################################################################################################################
# Combine Kriegstein and Costello
#############################################################################################################################################################

m1 <- MP_scores[, -c(1:4)]
m1 <- m1 %>%
  select(starts_with("MP_"), -MP_1_RP, -MP_11_MIC, -MP_12_LQ, -MP_16_GlioNeural)
m1 <- as.matrix(m1)
rownames(m1) <- MP_scores$CellID

m2 <- COSTELLO_scores[, -c(1:4)]
m2 <- as.matrix(m2)
rownames(m2) <- COSTELLO_scores$CellID

m3 <- KRIEG_scores[, -c(1:4)]
m3 <- as.matrix(m3)
rownames(m3) <- KRIEG_scores$CellID

m1 <- m1[rownames(m2), ]
m3 <- m3[rownames(m2), ]

m1 <- m1[, mps]
m2 <- m2[, costello]
m3 <- m3[, krieg]

colnames(m2) <- paste0("MATH_", colnames(m2))
colnames(m3) <- paste0("BHAD_", colnames(m3))

cor_m <- cor(m1, cbind(m2, m3))

x_labels <- c("MATH_Choroid" = "Choroid (Mathur)", "MATH_Astrocyte" = "Astrocyte (Mathur)",
              "MATH_DevelopmentalGSC" = "Developmental GSC (Mathur)", "MATH_EN.PFC1" = "EN PFC 1 (Mathur)",
              "MATH_EN.PFC2" = "EN PFC 2 (Mathur)", "MATH_InjuryResponseGSC" = "Injury Response GSC (Mathur)",
              "MATH_Interferon response" = "Interferon response (Mathur)", "MATH_IPC.div1" = "Divining IPC 1 (Mathur)",
              "MATH_IPC.div2" = "Divining IPC 2 (Mathur)", "MATH_MES2" = "MES2 (Mathur)",
              "MATH_Mesenchymal" = "Mesenchymal (Mathur)", "MATH_Neuron" = "Neuron (Mathur)",
              "MATH_OPC" = "OPC (Mathur)", "MATH_oRG" = "oRG (Mathur)", "MATH_RG.div1" = "Dividing RG 1 (Mathur)",
              "MATH_RG.div2" = "Dividing RG 2 (Mathur)",
              "BHAD_Glycolytic Progenitor_2" = "Glycolytic Progenitor 1 (Bhaduri)",
              "BHAD_Glycolytic Progenitor_34" = "Glycolytic Progenitor 2 (Bhaduri)",
              "BHAD_Protoplasmic Astrocyte_15" = "Protoplasmic Astrocyte (Bhaduri)",
              "BHAD_Immature Astrocyte_43" = "Immature Astrocyte 1 (Bhaduri)",
              "BHAD_Immature Astrocyte_41" = "Immature Astrocyte 2 (Bhaduri)",
              "BHAD_Radial Glia_38" = "Radial Glia 1 (Bhaduri)", "BHAD_Radial Glia_28" = "Radial Glia 2 (Bhaduri)",
              "BHAD_Radial Glia_19" = "Radial Glia 3 (Bhaduri)", "BHAD_Radial Glia_26" = "Radial Glia 4 (Bhaduri)",
              "BHAD_Radial Glia_27" = "Radial Glia 5 (Bhaduri)", "BHAD_Radial Glia_23" = "Radial Glia 6 (Bhaduri)",
              "BHAD_OPC_31" = "OPC 1 (Bhaduri)", "BHAD_OPC_11" = "OPC 2 (Bhaduri)",
              "BHAD_OPC_3" = "OPC 3 (Bhaduri)", "BHAD_OPC_39" = "OPC 4 (Bhaduri)",
              "BHAD_Mixed Progenitor/Neuron_6" = "Mixed Progenitor/Neuron 1 (Bhaduri)",
              "BHAD_Mixed Progenitor/Neuron_42" = "Mixed Progenitor/Neuron 2 (Bhaduri)",
              "BHAD_Mature IPC/Newborn Neuron_1" = "Mixed Progenitor/Neuron 3 (Bhaduri)",
              "BHAD_CGE iN_47" = "CGE iN (Bhaduri)", "BHAD_Neuron_45" = "Neuron (Bhaduri)",
              "BHAD_Dividing Neuron_14" = "Dividing Neuron 1 (Bhaduri)",
              "BHAD_Dividing Neuron_18" = "Dividing Neuron 2 (Bhaduri)",
              "BHAD_Dividing Neuron_30" = "Dividing Neuron 3 (Bhaduri)",
              "BHAD_Dividing Neuron_40" = "Dividing Neuron 4 (Bhaduri)",
              "BHAD_Dividing OPC_25" = "Dividing OPC (Bhaduri)",
              "BHAD_Dividing Progenitor_51" = "Dividing Progenitor 1 (Bhaduri)",
              "BHAD_Dividing Progenitor_52" = "Dividing Progenitor 2 (Bhaduri)")

y_labels <- c("MP_5_Hypoxia" = "Hypoxia", "MP_6_MES" = "MES-like", "MP_13_Cilia" = "Cilia-like",
              "MP_4_AC" = "AC-like", "MP_8_GPC" = "GPC-like", "MP_2_OPC" = "OPC-like",
              "MP_7_NPC" = "NPC-like", "MP_9_ExN" = "NEU-like", "MP_3_CC" = "Cell cycle")

colnames(cor_m) <- x_labels[colnames(cor_m)]
rownames(cor_m) <- y_labels[rownames(cor_m)]

hc1 <- fastcluster::hclust(d = dist(cor_m, "euclidean"), method = "complete")
hc2 <- fastcluster::hclust(d = dist(t(cor_m), "euclidean"), method = "complete")

dm <- melt(cor_m)

dm$Var1 <- factor(dm$Var1, hc1$labels[hc1$order])
dm$Var2 <- factor(dm$Var2, hc2$labels[hc2$order])

ggplot(dm, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Correlation", low = "dodgerblue", mid = "white", high = "red", limits = c(-.75, .75), oob = squish,
                       breaks = c(-.75, -.5, -.25, 0, .25, .5, .75), labels = c("-.75", "", "", "0", "", "", ".75")) +
  xlab("GBM states (other datasets)") +
  ylab("GBM states\n(this dataset)") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90),
                 axis.text.y = element_text(size = 16))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF4a - Fraction of cycling cells across states
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- mdata %>%
  filter(State %in% c("Cilia", "AC", "MES", "Hypoxia", "Stress", "GPC", "OPC", "NPC", "Neuron"))
table(d$State)
table(d$State, d$isCC)
table(d$State, d$isCC) / rowSums(table(d$State, d$isCC))

d$PvsR <- gsub(" ", "", d$PvsR)
d$PvsR[grep("Recurrent", d$PvsR)] <- "Recurrence"
table(d$PvsR)

d <- d %>%
  group_by(Sample, PvsR, State) %>%
  summarise(n = sum(isCC), N = n(), Freq = n / N, .groups = "drop") %>%
  filter(N >= 50)

d %>%
  group_by(State) %>%
  summarise(Mean = mean(Freq), SD = sd(Freq), n = n(), SE = SD / sqrt(n))

d <- d %>%
  mutate(State_fct = case_when(State == "Cilia" ~ "Cilia-like",
                               State == "AC" ~ "AC-like",
                               State == "MES" ~ "MES-like",
                               State == "Hypoxia" ~ "Hypoxia",
                               State == "Stress" ~ "Stress",
                               State == "GPC" ~ "GPC-like",
                               State == "OPC" ~ "OPC-like",
                               State == "NPC" ~ "NPC-like",
                               State == "Neuron" ~ "Neuron-like")) %>%
  mutate(State_fct = factor(State_fct, c("Cilia-like", "AC-like", "MES-like", "Hypoxia",
                                         "Stress", "GPC-like", "OPC-like", "NPC-like",
                                         "Neuron-like")))

d_ord <- d %>%
  group_by(State_fct) %>%
  summarise(Median = median(Freq), Mean = mean(Freq), .groups = "drop") %>%
  arrange(Median)

d$State_fct <- factor(as.character(d$State_fct), as.character(d_ord$State_fct))

ggboxplot(data = d, x = "State_fct", y = "Freq", add = "jitter", xlab = "State", ylab = "% cycling cells") +
  scale_y_continuous(labels = percent, limits = c(0, .4), oob = squish) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF4e - Lab contribution per MP
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

nmf_clusters <- malignant_nmf_metaprograms$clusters

res <- lapply(names(nmf_clusters), function(xname) {
  x <- nmf_clusters[[xname]]
  res <- sapply(strsplit(unique(x), split = "_"), function(y) y[[1]][1])
  table(sample_data$Lab[sample_data$Sample %in% res]) %>%
    melt() %>%
    mutate(cluster = xname)
})
res <- do.call(rbind, res) %>%
  as_tibble()

res$Var1 <- as.character(res$Var1)

n_samples_per_lab <- meta_data %>%
  select(Sample, Lab) %>%
  filter(!duplicated(Sample)) %>%
  group_by(Lab) %>%
  summarise(n = n(), .groups = "drop")

n_samples_per_lab <- setNames(n_samples_per_lab$n, n_samples_per_lab$Lab)

res$value <- res$value / n_samples_per_lab[res$Var1]

res$MP <- gsub("Cluster", "MP", res$cluster)

res <- res %>%
  filter(MP != "MP_16")

res$MP <- factor(res$MP, paste0("MP_", 1:15))

ggplot(res, aes(x = MP, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_discrete(name = "Lab") +
  scale_y_continuous(labels = percent) +
  xlab("Meta-program") +
  ylab("% contributing samples\n(proportional to lab-specific samples)") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 12)) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF4f - MC contribution per MP
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

nmf_clusters <- malignant_nmf_metaprograms$clusters

res <- lapply(names(nmf_clusters), function(xname) {
  x <- nmf_clusters[[xname]]
  res <- sapply(strsplit(unique(x), split = "_"), function(y) y[[1]][1])
  table(sample_data$Origin[sample_data$Sample %in% res]) %>%
    melt() %>%
    mutate(cluster = xname)
})
res <- do.call(rbind, res) %>%
  as_tibble()

res$Var1 <- as.character(res$Var1)

n_samples_per_origin <- meta_data %>%
  select(Sample, Origin) %>%
  filter(!duplicated(Sample)) %>%
  group_by(Origin) %>%
  summarise(n = n(), .groups = "drop")

n_samples_per_origin <- setNames(n_samples_per_origin$n, n_samples_per_origin$Origin)

res$value <- res$value / n_samples_per_origin[res$Var1]

res$MP <- gsub("Cluster", "MP", res$cluster)

res <- res %>%
  filter(MP != "MP_16")

res$MP <- factor(res$MP, paste0("MP_", 1:15))

ggplot(res, aes(x = MP, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_discrete(name = "Origin") +
  scale_y_continuous(labels = percent) +
  xlab("Meta-program") +
  ylab("% contributing samples\n(proportional to institution-specific samples)") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 12)) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF5a - Pathway MP heatmap
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

limits <- lengths(nmf_pathway_metaprograms$clusters) %>% cumsum()
limits <- c(1, limits[-length(limits)])

x1 <- limits[seq(from = 1, to = length(limits) - 1, by = 1)]
x2 <- limits[seq(from = 2, to = length(limits), by = 1)]
x1 <- factor(levels(nmf_pathway_metaprograms$nmf_intersect_meltI$Var1)[x1],
             levels(nmf_pathway_metaprograms$nmf_intersect_meltI$Var1))
x2 <- factor(levels(nmf_pathway_metaprograms$nmf_intersect_meltI$Var1)[x2],
             levels(nmf_pathway_metaprograms$nmf_intersect_meltI$Var1))
y1 <- x1; y2 <- x2
f <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)

ggplot(data = nmf_pathway_metaprograms$nmf_intersect_meltI,
       aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  # geom_tile() +
  geom_raster() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  geom_rect(data = f, aes(x = NULL,y = NULL, xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            color = "black", linetype = "dashed", fill = NA, size = .75, inherit.aes = F) +
  scale_y_discrete(limits = rev)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF5d - Hybrids complexity
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- mdata %>%
  mutate(isHybrid = mdata$HybridID != "Singular")

dm <- d %>%
  group_by(Sample, isHybrid) %>%
  summarise(comp = mean(nFeature_RNA), SD = sd(nFeature_RNA), N = n(), SE = SD / sqrt(N), .groups = "drop")

ggboxplot(data = dm, x = "isHybrid", y = "comp", add = "jitter") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

t.test(x = dm$comp[d$isHybrid], y = dm$comp[!d$isHybrid])

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF5e - Observed vs. expected hybrids frequency
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

obs_hybrid_freq <- mdata %>%
  group_by(HybridID) %>%
  summarise(n = n()) %>%
  mutate(N = sum(n), Freq = n / N)

exp_hybrid_freq <- mdata %>%
  group_by(State) %>%
  summarise(n = n()) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  arrange(State)

hybrid_ids <- unique(mdata$HybridID)

exp_hybrid <- lapply(hybrid_ids, function(hid) {
  hid1 <- strsplit(hid, split = "-")[[1]][1]
  hid2 <- strsplit(hid, split = "-")[[1]][2]
  res <- tibble(Hybrid = hid,
                ExpFreq = exp_hybrid_freq$Freq[exp_hybrid_freq$State == hid1] * exp_hybrid_freq$Freq[exp_hybrid_freq$State == hid2])
  return(res)
})
exp_hybrid <- do.call(rbind, exp_hybrid)

obs_hybrids <- melt(table(mdata$HybridID) / nrow(mdata)) %>%
  as_tibble()
colnames(obs_hybrids) <- c("Hybrid", "ObsFreq")

obs_vs_exp_hybrids <- exp_hybrid %>%
  left_join(obs_hybrids, by = "Hybrid")

obs_vs_exp_hybrids <- obs_vs_exp_hybrids[!is.na(obs_vs_exp_hybrids$ExpFreq) & !is.na(obs_vs_exp_hybrids$ObsFreq), ]

obs_vs_exp_hybrids$FC <- obs_vs_exp_hybrids$ObsFreq / obs_vs_exp_hybrids$ExpFreq
obs_vs_exp_hybrids$log2FC <- log2(obs_vs_exp_hybrids$FC)

obs_vs_exp_hybrids$Hybrid <- factor(obs_vs_exp_hybrids$Hybrid, obs_vs_exp_hybrids %>%
                                      arrange(log2FC) %>% pull(Hybrid))

technical_factor <- obs_vs_exp_hybrids %>%
  filter(log2FC < -1) %>%
  summarise(x = mean(FC)) %>%
  pull(x)

obs_vs_exp_hybrids$ExpFreqTechnical <- obs_vs_exp_hybrids$ExpFreq * technical_factor

obs_vs_exp_hybrids$log2FC_2 <- log2(obs_vs_exp_hybrids$ObsFreq / obs_vs_exp_hybrids$ExpFreqTechnical)

obs_vs_exp_hybrids$N <- obs_hybrid_freq$N[1]
obs_vs_exp_hybrids$NexpTechnical <- round(obs_vs_exp_hybrids$ExpFreqTechnical * obs_vs_exp_hybrids$N, 0)
obs_vs_exp_hybrids$Nexp <- round(obs_vs_exp_hybrids$ExpFreq * obs_vs_exp_hybrids$N, 0)
obs_vs_exp_hybrids$Nobs <- obs_vs_exp_hybrids$ObsFreq * obs_vs_exp_hybrids$N

obs_vs_exp_hybrids <- obs_vs_exp_hybrids %>%
  group_by(Hybrid) %>%
  mutate(pval = fisher.test(x = matrix(data = c(c(Nobs, N - Nobs),
                                                c(NexpTechnical, N - NexpTechnical)),
                                       byrow = T, nrow = 2), alternative = "greater")$p.value) %>%
  ungroup()

obs_vs_exp_hybrids$padj <- p.adjust(obs_vs_exp_hybrids$pval, "holm")

obs_vs_exp_hybrids$sigRes <- obs_vs_exp_hybrids$padj < .05

obs_vs_exp_hybrids$Nobs_per1000 <- round(obs_vs_exp_hybrids$ObsFreq * 1000, 0)
obs_vs_exp_hybrids$Nexp_per1000 <- round(obs_vs_exp_hybrids$ExpFreq * 1000, 0)
obs_vs_exp_hybrids$NexpTechnical_per1000 <- round(obs_vs_exp_hybrids$ExpFreqTechnical * 1000, 0)

obs_vs_exp_hybrids <- obs_vs_exp_hybrids %>%
  group_by(Hybrid) %>%
  mutate(pval = fisher.test(x = matrix(data = c(c(Nobs_per1000, 1000 - Nobs_per1000),
                                                c(NexpTechnical_per1000, 1000 - NexpTechnical_per1000)),
                                       byrow = T, nrow = 2), alternative = "greater")$p.value) %>%
  ungroup()

obs_vs_exp_hybrids$padj <- p.adjust(obs_vs_exp_hybrids$pval, "fdr")

obs_vs_exp_hybrids$sigRes <- obs_vs_exp_hybrids$padj < .05

obs_vs_exp_hybrids$Label <- ifelse(obs_vs_exp_hybrids$sigRes, "*", "")

ggplot(obs_vs_exp_hybrids, aes(x = Hybrid, y = log2FC_2, fill = sigRes)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = Label, y = 6), color = "black", size = 10) +
  scale_fill_brewer(name = "", palette = "Set1") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 1) +
  xlab("Hybrid pair") +
  ylab("Observed vs. expected\ntechnical hgybrids [log2]") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90), legend.position = "none") +
  theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 3 - State-controlled baseline profiles
#
# This section generates the state-controlled baseline profiles that characterize the inter-tumor heterogenity based on genes that
# recur across state-specific pseudo-bulk profiles
#
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
# Generate the PB profiles
####################################################################################################################################

d <- mdata
length(table(d$Sample))

# Generate profiles from the 7 principal cellular states
d_state_stats <- d %>%
  filter(State %in% c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")) %>%
  group_by(Patient, Timepoint, Sample, State) %>%
  summarise(n = n())

n_cells <- 25

states <- c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")

d_state_stats <- d_state_stats %>%
  filter(n >= n_cells)

# Remove "junk genes" (i.e. pseudo-genes, antisense genes etc)
junk_genes <- c(rownames(umi_data_all[[1]])[grep("\\.", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("-AS*", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("LINC", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^RP[S|L]", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^MT-", rownames(umi_data_all[[1]]))])

valid_genes <- rownames(umi_data_all[[1]])[rownames(umi_data_all[[1]]) %ni% junk_genes]

# Generate the pseudo-bulk profiles
state_pb_profiles <- lapply(unique(d_state_stats$State), function(st) {
  
  print(st)
  
  d_state_stats_tmp <- d_state_stats %>%
    filter(State == st)
  
  res <- lapply(d_state_stats_tmp$Sample, function(sname) {
    
    m <- umi_data_all[[sname]]
    
    d_tmp <- d %>%
      filter(Sample == sname, State == st)
    
    m <- m[valid_genes, d_tmp$CellID]
    
    m <- umi2upm(m)
    
    m <- log2(rowMeans(m) + 1)
    
    return(m)
  })
  res <- do.call(cbind, res)
  colnames(res) <- paste0(d_state_stats_tmp$Sample, "_", st)
  
  return(res)
})
names(state_pb_profiles) <- unique(d_state_stats$State)

####################################################################################################################################
# Select genes for the analysis
####################################################################################################################################

complexity <- apply(do.call(cbind, state_pb_profiles), 2, function(x) sum(x > 0))
n_cells <- setNames(d_state_stats$n, paste0(d_state_stats$Sample, "_", d_state_stats$State))

# Remove genes that don't have an entrezid
valid_genes <- tibble(Gene = rownames(state_pb_profiles[[1]]),
                      Entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, Gene, 'ENTREZID', 'SYMBOL'))
valid_genes <- valid_genes %>%
  filter(!is.na(Entrez))

rm_list <- lapply(state_pb_profiles, rowMeans)
rm_list <- lapply(rm_list, function(x) x[names(x) %in% valid_genes$Gene])

rm_list <- lapply(rm_list, function(x) x[x > 4])

rm <- rowMeans(do.call(cbind, state_pb_profiles))
rm <- rm[names(rm) %in% valid_genes$Gene]

# Variance is computed within each state to capture the genes that change within state (rather than across states)
rv <- lapply(state_pb_profiles, function(x) apply(x, 1, var))
rv <- do.call(cbind, rv)
rv <- apply(rv, 1, median)
rv <- rv[names(rv) %in% valid_genes$Gene]
rv <- rv[names(rm)]

sg <- tibble(Gene = names(rv), V = rv, M = rm)

# Select genes with mean expression > 1 abd variance > 2.5 (across all pseudo-bulk profiles)
sg <- sg %>%
  filter(M > 1, V > 2.5)
dim(sg)

sg <- sg$Gene

####################################################################################################################################
# Compute PCA per state to detect the genes contributing most variance within each state
####################################################################################################################################

pca_res_list <- lapply(names(state_pb_profiles), function(st) {
  
  print(st)
  
  pbm <- state_pb_profiles[[st]][sg, ]
  dim(pbm)
  
  # Center within state  
  pbm <- rowcenter(pbm)
  
  pca_res <- pcaMethods::pca(object = t(pbm), nPcs = 10, scale = "none", center = F)
  
  return(pca_res)
})
names(pca_res_list) <- names(state_pb_profiles)

pca_scores_tbl <- do.call(rbind,
                          lapply(pca_res_list, function(x) as_tibble(pcaMethods::scores(x), rownames = "ID")))
pca_scores_tbl$Sample <- sapply(strsplit(pca_scores_tbl$ID, "_"), function(x) x[[1]][1])
pca_scores_tbl$State <- sapply(strsplit(pca_scores_tbl$ID, "_"), function(x) x[[2]][1])

####################################################################################################################################
# Generate two signatures from each (State, PC)
####################################################################################################################################

npcs <- 3
n_genes <- 50

pca_loadings_tbl <- lapply(names(pca_res_list), function(st) {
  
  print(st)
  
  pca_res <- pca_res_list[[st]]
  
  lds <- pcaMethods::loadings(pca_res)
  dim(lds)
  
  lds_res <- lapply(1:npcs, function(i) {
    res <- tibble(V1 = lds[, i] %>% sort(decreasing = T) %>% head(n_genes) %>% names(),
                  V2 = lds[, i] %>% sort(decreasing = F) %>% head(n_genes) %>% names())
    
    colnames(res) <- paste0(st, "_PC", i, c("_HIGH", "_LOW"))
    
    return(res)
  })
  lds_res <- do.call(cbind, lds_res)
  
  return(lds_res)
})
pca_loadings_tbl <- do.call(cbind, pca_loadings_tbl)

####################################################################################################################################
# Cluster the signatures using Jaccard index as similarity metric
####################################################################################################################################

jac_m <- jaccard(pca_loadings_tbl)

pca_ccp <- ConsensusClusterPlus::ConsensusClusterPlus(d = dist(jac_m, "euclidean"), maxK = 10, reps = 1000, innerLinkage = "average", finalLinkage = "average", plot = "pdf")
pca_k <- 5
clusters <- as.factor(pca_ccp[[pca_k]]$consensusClass)

hc <- pca_ccp[[pca_k]]$consensusTree

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF6a - Jaccard similiarity matrix of PCA signatures
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

Heatmap(jac_m, use_raster = F, raster_by_magick = T,
        col = circlize::colorRamp2(c(0, .03, .06, .15, .2, .25), c("white", "white", "gold", "orange", "red", "darkred")),
        show_row_names = T, show_column_names = T,
        cluster_columns = hc, cluster_rows = hc,
        clustering_distance_rows = "euclidean", clustering_method_rows = "average",
        clustering_distance_columns = "euclidean", clustering_method_columns = "average",
        column_title = "Programs", column_title_gp = gpar(fontsize = 20),
        row_title = "Programs", row_title_gp = gpar(fontsize = 20), 
        row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(fontsize = 14),
        top_annotation = HeatmapAnnotation(Cluster = clusters,
                                           col = list(Cluster = setNames(brewer_pal(palette = "Spectral")(length(levels(clusters))), levels(clusters))),
                                           show_annotation_name = FALSE, #show_legend = FALSE,
                                           annotation_legend_param = list(title_gp = gpar(fontsize = 20),
                                                                          labels_gp = gpar(fontsize = 14), border = "black",
                                                                          grid_height = unit(7, "mm"))),
        right_annotation = rowAnnotation(Cluster = clusters,
                                         col = list(Cluster = setNames(brewer_pal(palette = "Spectral")(length(levels(clusters))), levels(clusters))),
                                         show_annotation_name = FALSE, show_legend = FALSE,
                                         annotation_legend_param = list(title_gp = gpar(fontsize = 20),
                                                                        labels_gp = gpar(fontsize = 14), border = "black",
                                                                        grid_height = unit(7, "mm"))),
        heatmap_legend_param = list(at = c(0, .05, .1, .15, .2, .25), title = "Jaccard",
                                    labels = c("0%", "5%", "10%", "15%", "20%", "25%"),
                                    title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 14), title_position = "topleft", border = "black",
                                    grid_height = unit(7, "mm")),
        column_title_side = "bottom",
        row_gap = unit(0,"mm"), cluster_column_slices = T, cluster_row_slices = F)

####################################################################################################################################
# Generate consensus signatures from the 5 clusters
####################################################################################################################################

pca_gene_sigs_ranked <- lapply(levels(clusters), function(cl) {
  n <- table(clusters)[cl]
  min_n <- round(max(.25*n, 3))
  print(min_n)
  gs <- table(unlist(pca_loadings_tbl[, pca_ccp[[pca_k]]$consensusClass == cl])) %>% sort()
  gs <- gs[gs >= min_n]
  return(gs)
})
names(pca_gene_sigs_ranked) <- paste0("C", levels(clusters))

pca_gene_sigs <- lapply(pca_gene_sigs_ranked, names)

####################################################################################################################################
# Gene-set enrichment to aid in annotating the consensus signatures
####################################################################################################################################

ego_mps <- lapply(pca_gene_sigs, function(x) {
  clusterProfiler::enrichGO(gene          = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL'),
                            universe      = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, names(rm), 'ENTREZID', 'SYMBOL'),
                            OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            minGSSize     = 50,
                            maxGSSize     = 300,
                            pool          = T,
                            readable      = TRUE)
})

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF6b - Go enrichment interaction graphs
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

edo_mp <- enrichplot::pairwise_termsim(ego_mps[[1]])
enrichplot::emapplot(edo_mp, showCategory = 10)

edo_mp <- enrichplot::pairwise_termsim(ego_mps[[3]])
enrichplot::emapplot(edo_mp, showCategory = 10)

edo_mp <- enrichplot::pairwise_termsim(ego_mps[[5]])
enrichplot::emapplot(edo_mp, showCategory = 10)

####################################################################################################################################
# Create a Gene x PB matrix containing all consensus genes to facilitate scoring
####################################################################################################################################

pca_gene <- unlist(pca_gene_sigs, recursive = F)  %>% unique() %>% unname()

pbm_list <- lapply(names(state_pb_profiles), function(st) {
  
  print(st)
  
  pbm <- state_pb_profiles[[st]][pca_gene, ]
  dim(pbm)
  
  pbm <- rowcenter(pbm)
  
  return(pbm)
})
pbm_list <- do.call(cbind, pbm_list)

####################################################################################################################################
# Score the PB profiles for the consensus signatures and create the feature_scores table
####################################################################################################################################

feature_scores <- tibble(ID = colnames(pbm_list),
                         C1 = colMeans(pbm_list[pca_gene_sigs$C1, ]),
                         C2 = colMeans(pbm_list[pca_gene_sigs$C2, ]),
                         C3 = colMeans(pbm_list[pca_gene_sigs$C3, ]),
                         C4 = colMeans(pbm_list[pca_gene_sigs$C4, ]),
                         C5 = colMeans(pbm_list[pca_gene_sigs$C5, ]),
                         Diff = C3 - C1)
feature_scores$Sample <- sapply(strsplit(feature_scores$ID, split = "_"), function(x) x[[1]][1])
feature_scores$State <- sapply(strsplit(feature_scores$ID, split = "_"), function(x) x[[2]][1])

feature_scores <- feature_scores %>%
  left_join(sample_data %>%
              select(Sample, Patient, Timepoint),
            by = "Sample")

feature_scores$Complexity <- complexity[feature_scores$ID]
feature_scores$Ncell <- n_cells[feature_scores$ID]

feature_scores$SID <- paste0(feature_scores$Patient, feature_scores$Timepoint)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 3b - Inter vs. intra tumor heterogenity
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

mp_list <- c(MP_list_named,
             list(BP_ECM = pca_gene_sigs$C1,
                  BP_Neuronal = pca_gene_sigs$C3,
                  BP_Glial = pca_gene_sigs$C5))

mp_vs_bp_scores <- lapply(names(mp_list), function(mpname) {
  print(mpname)
  mp <- mp_list[[mpname]]
  res <- lapply(state_pb_profiles, function(pb) {
    pb <- pb[rownames(pb) %in% mp, ]
    tibble(ID = colnames(pb),
           Sample = sapply(strsplit(ID, split = "_"), function(x) x[[1]][1]),
           State = sapply(strsplit(ID, split = "_"), function(x) x[[2]]),
           MP = mpname, Score = colMeans(pb))
  })
  res <- do.call(rbind, res)
  return(res)
})
mp_vs_bp_scores <- do.call(rbind, mp_vs_bp_scores)

x <- mp_vs_bp_scores %>%
  group_by(Sample, MP) %>%
  summarise(SD = sd(Score), .groups = "drop") %>%
  filter(!is.na(SD))

x <- x %>%
  group_by(MP) %>%
  summarise(Intra = mean(SD), .groups = "drop")

y <- mp_vs_bp_scores %>%
  group_by(State, MP) %>%
  summarise(SD = sd(Score), .groups = "drop")

y <- y %>%
  group_by(MP) %>%
  summarise(Inter = mean(SD), .groups = "drop")

d <- x %>%
  left_join(y, by = "MP")

d <- d %>%
  filter(MP %in% c("MP_2_OPC", "MP_4_AC", "MP_5_Hypoxia", "MP_6_MES", "MP_7_NPC", "MP_8_GPC", "MP_9_ExN",
                   "BP_ECM", "BP_Neuronal", "BP_Glial"))

d$Type <- sapply(strsplit(d$MP, split = "_"), function(x) x[[1]][1])

ggplot(d, aes(x = Intra, y = Inter)) +
  geom_point(aes(fill = Type), color = "black", shape = 21, size = 8) +
  geom_text(aes(label = MP), size = 8) +
  scale_fill_brewer(name = "", palette = "Set1", labels = c("BP" = "Baseline profile", "MP" = "Meta-program")) +
  scale_x_continuous(limits = c(.5, 1)) +
  scale_y_continuous(limits = c(.5, 1.3)) +
  xlab("Intra-tumor variability") +
  ylab("Inter-tumor variability") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

t.test(d$Inter[d$Type == "BP"], d$Inter[d$Type == "MP"], var.equal = T, alternative = "greater")
t.test(d$Intra[d$Type == "BP"], d$Intra[d$Type == "MP"], var.equal = T, alternative = "less")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 3c - BP-ECM - BP-Neuronal difference for each sample
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

feature_scores$Ord <- factor(feature_scores$ID,
                             feature_scores %>%
                               arrange(Diff) %>%
                               pull(ID))

gene_ord <- cor(t(pbm_list[pca_gene, ]), feature_scores$Diff) %>%
  as_tibble(rownames = "Gene") %>%
  arrange(V1)

sample_ord <- feature_scores %>%
  group_by(SID) %>%
  summarise(Mean = mean(Diff)) %>%
  arrange(Mean)

feature_scores$Sample_ord <- factor(feature_scores$SID, sample_ord$SID)

ggplot(feature_scores, aes(x = Sample_ord, y = Diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = State), color = "black", shape = 21, size = 3, width = .2, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Samples") +
  ylab("Score difference (BP-ECM - BP-Neuronal)") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF6c - Within-sample variability for each of the scores
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ggplot(feature_scores %>%
         mutate(Sample_ord = factor(Sample, feature_scores %>%
                                      group_by(Sample) %>%
                                      summarise(Mean = mean(C1)) %>%
                                      arrange(Mean) %>%
                                      pull(Sample))),
       aes(x = Sample_ord, y = C1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = State), color = "black", shape = 21, size = 1, width = .2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("BP-ECM score") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_line(), legend.position = "none") +
  ggplot(feature_scores %>%
           mutate(Sample_ord = factor(Sample, feature_scores %>%
                                        group_by(Sample) %>%
                                        summarise(Mean = mean(C3)) %>%
                                        arrange(Mean) %>%
                                        pull(Sample))),
         aes(x = Sample_ord, y = C3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = State), color = "black", shape = 21, size = 1, width = .2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("BP-Neuronal score") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_line(), legend.position = "none") +
  ggplot(feature_scores %>%
           mutate(Sample_ord = factor(Sample, feature_scores %>%
                                        group_by(Sample) %>%
                                        summarise(Mean = mean(C5)) %>%
                                        arrange(Mean) %>%
                                        pull(Sample))),
         aes(x = Sample_ord, y = C5)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = State), color = "black", shape = 21, size = 1, width = .2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Samples") +
  ylab("BP-Glial score") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_line(), legend.position = "none") +
  plot_layout(ncol = 1)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF6d - Compute and plot the observed vs. expected within-sample variability in BP-Neuronal - BP-ECM difference
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

obs_sample_sd <- feature_scores %>%
  group_by(Sample) %>%
  summarise(n = n(), SD = sd(Diff)) %>%
  filter(n >= 2) %>%
  ungroup()

exp_sample_sd <- lapply(1:100, function(i) {
  res <- feature_scores %>%
    sample_n(3) %>%
    summarise(n = n(), SD = sd(Diff))
  return(res)
})
exp_sample_sd <- do.call(rbind, exp_sample_sd)

obs_vs_exp_sample_sd <- rbind(obs_sample_sd %>% mutate(DataType = paste0("Observed")) %>% select(SD, DataType),
                              exp_sample_sd %>% select(SD) %>% mutate(DataType = "Expected"))

ggplot(obs_vs_exp_sample_sd, aes(x = DataType, y = SD)) + 
  ggdist::stat_halfeye(aes(fill = DataType), color = "black", adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA) + 
  geom_boxplot(width = .12, outlier.color = NA) +
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .3) +
  coord_cartesian(xlim = c(1.2, NA)) +
  scale_fill_brewer(name = "", palette = "Set1") +
  xlab("") +
  ylab("BPP-Neuronal - BP-ECM\ndifference intra-sample SD") +
  scale_y_continuous(limits = c(0, 3), oob = squish) +
  theme_gbm_pvsr() +
  theme(panel.grid.major = element_line())

####################################################################################################################################
# 3-axis coordinates for C1/3/5 scores to depict inter-tumor heterogeneity
####################################################################################################################################

feature_scores$Lineage <- apply(feature_scores[, c("C1", "C3")], 1, max)
feature_scores$Stemness <- apply(feature_scores[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

feature_scores$LineagePlot <- apply(feature_scores[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.25)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF6e - 3-axis per state
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ggplot(feature_scores, aes(x = LineagePlot, y = C5)) +
  geom_point(aes(fill = State), color = "black", shape = 22, size = 3, alpha = .75) +
  geom_smooth(method = "loess", color = "black", se = F, size = 2, span = .75) +
  scale_fill_brewer(palette = "Set3") +
  xlab("ECM <-> Neuronal") +
  ylab("-> Glial") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
# Generate feature scores per sample the and 3-axis coordinates
####################################################################################################################################

feature_scores_per_sample <- feature_scores %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(C1 = mean(C1), C2 = mean(C2), C3 = mean(C3), C4 = mean(C4), C5 = mean(C5)) %>%
  ungroup()

feature_scores_per_sample$Lineage <- apply(feature_scores_per_sample[, c("C1", "C3")], 1, max)
feature_scores_per_sample$Stemness <- apply(feature_scores_per_sample[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

feature_scores_per_sample$LineagePlot <- apply(feature_scores_per_sample[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.25)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

####################################################################################################################################
# Generate a color code for each axis to facilitate using 3 color scales
####################################################################################################################################

feature_scores_per_sample$MaxAxis <- apply(feature_scores_per_sample[, c("C1", "C3", "C5")], 1,  function(x) c("C1", "C3", "C5")[which.max(x)])

feature_scores_per_sample$C1_rescaled <- rescale(feature_scores_per_sample$C1, c(0, .3))
feature_scores_per_sample$C3_rescaled <- rescale(feature_scores_per_sample$C3, c(.35, .65))
feature_scores_per_sample$C5_rescaled <- rescale(feature_scores_per_sample$C5, c(.7, 1))

feature_scores_per_sample <- feature_scores_per_sample %>%
  mutate(axisColorCode = case_when(MaxAxis == "C1" ~ C1_rescaled,
                                   MaxAxis == "C3" ~ C3_rescaled,
                                   MaxAxis == "C5" ~ C5_rescaled,
                                   TRUE ~ NA))

axis_color_vec <- c(RColorBrewer::brewer.pal(n = 9, name = "Greens"),
                    RColorBrewer::brewer.pal(n = 9, name = "Purples"),
                    RColorBrewer::brewer.pal(n = 9, name = "Reds"))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 3d (upper panel) - 3-axis per sample
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

feature_scores_per_sample$ID <- paste0(feature_scores_per_sample$Patient, feature_scores_per_sample$Timepoint)

ggplot(feature_scores_per_sample %>%
         mutate(Label = ifelse(ID %in% c("P64T2", "P66T2", "P51T1", "P47T1", "P76T2", "P37T2"), ID, "")),
       aes(x = LineagePlot, y = C5)) +
  geom_point(aes(fill = axisColorCode), color = "black", shape = 22, size = 5) +
  geom_text(aes(label = Label), size = 8, hjust = "left", vjust = "outward") +
  geom_smooth(method = "loess", color = "black", size = 1, se = F) +
  scale_fill_gradientn(name = "", colours = axis_color_vec, guide = "colorsteps", rescaler = function(x, ...) { return(x) },
                       breaks = seq(0, 1, .05), labels = NULL,
                       values = c(0, .1, .125, .15, .175, .2, .225, .25, .3,
                                  .35, .45, .475, .5, .525, .55, .575, .6, .65,
                                  .7, .8, .825, .85, .875, .9, .925, .95, 1)) +
  scale_x_continuous(limits = c(-3, 3), oob = squish, breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(-3, 2), oob = squish) +
  xlab("ECM <-> Neuronal") +
  ylab("-> Glial") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
# Score the cells for the 3 consensus signatures
####################################################################################################################################

m <- umi_data[, mdata$CellID]

m <- umi2upm(m)

rm <- log2(rowMeans(m) + 1)
sg <- names(rm[rm > 4])

m <- m[sg, ]

scp_scores <- lapply(unique(mdata$State), function(st) {
  
  print(st)
  
  d_st <- mdata %>%
    filter(State == st)
  
  print(nrow(d_st))
  
  m_st <- m[, d_st$CellID]
  m_st <- log2(m_st / 10 + 1)
  m_st <- as.matrix(m_st)
  
  sigScores(m = m_st, sigs = pca_gene_sigs) %>%
    as_tibble(rownames = "CellID")
})
scp_scores <- do.call(rbind, scp_scores)

scp_scores$Lineage <- apply(scp_scores[, c("C1", "C3")], 1, max)
scp_scores$Stemness <- apply(scp_scores[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

scp_scores$LineagePlot <- apply(scp_scores[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.5)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

scp_scores <- scp_scores %>%
  left_join(mdata %>%
              select(State, CellID), by = "CellID")

####################################################################################################################################
# Color-code the cells to facilitate plotting the butterflies
####################################################################################################################################

scp_scores$MaxAxis <- apply(scp_scores[, c("C1", "C3", "C5")], 1,  function(x) c("C1", "C3", "C5")[which.max(x)])

scp_scores$C1_rescaled <- rescale(scp_scores$C1, c(0, .3))
scp_scores$C3_rescaled <- rescale(scp_scores$C3, c(.35, .65))
scp_scores$C5_rescaled <- rescale(scp_scores$C5, c(.7, 1))

scp_scores <- scp_scores %>%
  mutate(axisColorCode = case_when(MaxAxis == "C1" ~ C1_rescaled,
                                   MaxAxis == "C3" ~ C3_rescaled,
                                   MaxAxis == "C5" ~ C5_rescaled,
                                   TRUE ~ NA))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 3d (lower panel) - 3-axis per sample
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- mdata

d <- d %>%
  ungroup() %>%
  left_join(scp_scores %>%
              select(CellID, MaxAxis, axisColorCode),
            by = "CellID")

ggplot(data = d %>%
         filter(ID %in% c("P64T2", "P66T2", "P51T1", "P47T1", "P76T2", "P37T2")) %>%
         mutate(ID_factor = factor(ID, c("P64T2", "P66T2", "P51T1", "P47T1", "P76T2", "P37T2"))),
       mapping = aes(x = Dx, y = Dy)) +
  facet_grid(cols = vars(ID_factor)) +
  geom_point(aes(fill = axisColorCode), color = "black", shape = 21, size = 3) +
  scale_fill_gradientn(name = "", colours = axis_color_vec, guide = "colorsteps", rescaler = function(x, ...) { return(x) },
                       breaks = seq(0, 1, .05), labels = NULL,
                       values = c(0, .1, .125, .15, .175, .2, .225, .25, .3,
                                  .35, .475, .5, .525, .55, .575, .6, .65,
                                  .7, .8, .825, .85, .875, .9, .925, .95, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, 1.5)) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(legend.position = "none")

####################################################################################################################################
# Generate the SCP scores for all samples (include also samples with few malignant cells)
####################################################################################################################################

d <- mdata

# Generate profiles from the 7 principal cellular states
d_state_stats <- d %>%
  filter(State %in% c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")) %>%
  group_by(Patient, Timepoint, Sample, State) %>%
  summarise(n = n())

states <- c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")

d_state_stats <- d_state_stats %>%
  filter(n >= 5)

# Remove "junk genes" (i.e. pseudo-genes, antisense genes etc)
junk_genes <- c(rownames(umi_data_all[[1]])[grep("\\.", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("-AS*", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("LINC", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^RP[S|L]", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^MT-", rownames(umi_data_all[[1]]))])

valid_genes <- rownames(umi_data_all[[1]])[rownames(umi_data_all[[1]]) %ni% junk_genes]

# Generate the pseudo-bulk profiles
state_pb_profiles_all <- lapply(unique(d_state_stats$State), function(st) {
  
  print(st)
  
  d_state_stats_tmp <- d_state_stats %>%
    filter(State == st)
  
  res <- lapply(d_state_stats_tmp$Sample, function(sname) {
    
    m <- umi_data_all[[sname]]
    
    d_tmp <- d %>%
      filter(Sample == sname, State == st)
    
    m <- m[valid_genes, d_tmp$CellID]
    
    m <- umi2upm(m)
    
    m <- log2(rowMeans(m) + 1)
    
    return(m)
  })
  res <- do.call(cbind, res)
  colnames(res) <- paste0(d_state_stats_tmp$Sample, "_", st)
  
  return(res)
})
names(state_pb_profiles_all) <- unique(d_state_stats$State)

pbm_all_list <- lapply(names(state_pb_profiles_all), function(st) {
  
  print(st)
  
  pbm <- state_pb_profiles_all[[st]][pca_gene, ]
  dim(pbm)
  
  pbm <- rowcenter(pbm)
  
  return(pbm)
})
pbm_all_list <- do.call(cbind, pbm_all_list)

feature_scores_all <- tibble(ID = colnames(pbm_all_list),
                             C1 = colMeans(pbm_all_list[pca_gene_sigs$C1, ]),
                             C3 = colMeans(pbm_all_list[pca_gene_sigs$C3, ]),
                             C5 = colMeans(pbm_all_list[pca_gene_sigs$C5, ]))
feature_scores_all$Sample <- sapply(strsplit(feature_scores_all$ID, split = "_"), function(x) x[[1]][1])
feature_scores_all$State <- sapply(strsplit(feature_scores_all$ID, split = "_"), function(x) x[[2]][1])

feature_scores_all <- feature_scores_all %>%
  left_join(sample_data %>%
              select(Sample, Patient, Timepoint),
            by = "Sample")

feature_scores_all_per_sample <- feature_scores_all %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(C1 = mean(C1), C3 = mean(C3), C5 = mean(C5)) %>%
  ungroup()

feature_scores_all_per_sample$Lineage <- apply(feature_scores_all_per_sample[, c("C1", "C3")], 1, max)
feature_scores_all_per_sample$Stemness <- apply(feature_scores_all_per_sample[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

feature_scores_all_per_sample$LineagePlot <- apply(feature_scores_all_per_sample[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.25)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

feature_scores_all_per_sample <- feature_scores_all_per_sample %>%
  group_by(Sample) %>%
  mutate(SCP = c("SCP-ECM", "SCP-Neuronal", "SCP-Glial")[which.max(c(C1, C3, C5))]) %>%
  ungroup()

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF7b - Correlation between snRNAseq and Spatial BPs
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

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

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 5 - GBM ecosystems
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 5b - Composition clusters
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- meta_data

d$CellType <- d$CellType_reduced

d_stats <- d %>%
  group_by(Sample, Patient, Timepoint, CellType, ID) %>%
  summarise(n = n()) %>%
  group_by(Sample, Patient, Timepoint, ID) %>%
  mutate(N = sum(n), Freq = n / N)

d_stats_m <- acast(d_stats, formula = ID ~ CellType, value.var = "Freq")
d_stats_m[is.na(d_stats_m)] <- 0

# Perform hierarchical clustering to aid in separating samples
hc <- fastcluster::hclust(d = dist(d_stats_m[, c("Malignant", "TAM", "Oligodendrocyte",
                                                 "Endothel", "Pericyte", "Astrocyte", "Excitatory neuron")], "euclidean"), method = "average")

d_stats$SampleCluster <- cutree(hc, h = .25)[as.character(d_stats$ID)] %>% as.character()

d_stats$ID <- factor(as.character(d_stats$ID), hc$labels[hc$order])

d_stats_m <- acast(data = d_stats, formula = ID ~ CellType, value.var = "Freq")
d_stats_m[is.na(d_stats_m)] <- 0

d_stats_m <- melt(d_stats_m) %>%
  as_tibble()
colnames(d_stats_m) <- c("ID", "CellType", "Freq")

d_stats_m$SampleCluster <- cutree(hc, h = .25)[as.character(d_stats_m$ID)]

cluster_type <- d_stats_m %>%
  group_by(ID) %>%
  summarise(Purity = Freq[CellType == "Malignant"], OC = Freq[CellType == "Oligodendrocyte"], TAM = Freq[CellType == "TAM"],
            GN = Freq[CellType == "Oligodendrocyte"] + Freq[CellType == "Astrocyte"] + Freq[CellType == "Excitatory neuron"] + Freq[CellType == "Inhibitory neuron"],
            EN = Freq[CellType == "Endothel"] + Freq[CellType == "Pericyte"])

cluster_type <- cluster_type %>%
  mutate(ClusterType = case_when(Purity > .75 ~"HP",
                                 Purity > .5 ~ "IP",
                                 TAM > .4 ~ "LP - TAM",
                                 OC > .4 ~ "LP - OC",
                                 GN > .4 ~ "LP - GN",
                                 # EN > .3 ~ "LP - EC",
                                 TRUE ~ "Mixed"))

d_stats$SampleCluster_int <- d_stats$SampleCluster

sc_vec <- setNames(cluster_type$ClusterType, cluster_type$ID)

d_stats$SampleCluster <- sc_vec[as.character(d_stats$ID)]
d_stats$SampleCluster <- factor(d_stats$SampleCluster, c("HP", "IP", "LP - TAM", "LP - OC", "LP - GN", "Mixed"))

d_stats %>%
  ggplot(aes(x = ID, y = Freq, fill = CellType)) +
  facet_grid(cols = vars(SampleCluster), scales = "free", space = "free_x") +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(name = "Cell type", values = celltype_color_vec_reduced) +
  xlab("") +
  ylab("Proportion") +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(panel.spacing = unit(1, "lines")) +
  theme(axis.text.x = element_text(size = 12, angle = 90))

comp_cluster_stats <- d_stats

comp_cluster_data <- comp_cluster_stats %>%
  group_by(Sample, ID, Patient, Timepoint) %>%
  summarise(CompCluster = first(SampleCluster), CompCluster_int = first(SampleCluster_int),
            SI_Comp = shannon_index(Freq),
            .groups = "drop")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 5c - Malignant state clusters
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

valid_samples <- mdata %>%
  group_by(Sample) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n >= 50)

d_m <- mdata %>%
  filter(Sample %in% valid_samples$Sample)

d_m_stats <- d_m %>%
  group_by(Sample, Patient, Timepoint, State, ID) %>%
  summarise(n = n()) %>%
  group_by(Sample, Patient, Timepoint, ID) %>%
  mutate(N = sum(n), Freq = n / N)

d_m_stats <- rbind(d_m_stats,
                   d_m %>%
                     group_by(Sample, Patient, Timepoint, isCC, ID) %>%
                     summarise(n = n()) %>%
                     group_by(Sample, Patient, Timepoint, ID) %>%
                     mutate(N = sum(n), Freq = n / N) %>%
                     filter(isCC == T) %>%
                     mutate(isCC = "Cycling") %>%
                     rename(State = isCC))

d_m_stats_m <- acast(data = d_m_stats, formula = ID ~ State, value.var = "Freq")
d_m_stats_m[is.na(d_m_stats_m)] <- 0
d_m_stats <- melt(d_m_stats_m) %>%
  as_tibble()
colnames(d_m_stats) <- c("ID", "State", "Freq")

d_m_stats <- d_m_stats %>%
  left_join(d_m %>%
              group_by(Sample, Patient, Timepoint, ID) %>%
              summarise(n = n()) %>%
              ungroup() %>%
              select(-n),
            by = "ID")

d_m_stats$ID <- factor(as.character(d_m_stats$ID), hc$labels[hc$order])

d_m_stats$State <- factor(as.character(d_m_stats$State),
                          c("Cilia", "AC", "MES", "Hypoxia", "Stress", "GPC", "OPC", "NPC", "Neuron", "Unresolved", "Cycling") %>% rev())

hc <- fastcluster::hclust(d = dist(d_m_stats_m, "euclidean"), method = "average")

d_m_stats$SampleCluster <- cutree(hc, k = 14)[as.character(d_m_stats$ID)] %>% as.character()

d_m_stats$ID <- factor(as.character(d_m_stats$ID), hc$labels[hc$order])

d_m_stats %>%
  ggboxplot(x = "State", y = "Freq", facet.by = "SampleCluster", add = "jitter") +
  theme(axis.text.x = element_text(angle = 90))

cluster_type <- d_m_stats %>%
  mutate(ID = as.character(ID)) %>%
  group_by(ID) %>%
  mutate(ClusterType = case_when(max(Freq) > .25 ~ as.character(State[which.max(Freq)]),
                                 TRUE ~ "Mixed"))
cluster_type %>%
  ggboxplot(x = "State", y = "Freq", facet.by = "ClusterType", add = "jitter") +
  theme(axis.text.x = element_text(angle = 90))

cluster_type <- cluster_type %>%
  mutate(ClusterType = case_when(ClusterType %in% c("MES", "Hypoxia") ~ "MES/Hypoxia",
                                 ClusterType %in% c("OPC", "NPC") ~ "OPC/NPC",
                                 ClusterType %in% c("Mixed", "Unresolved") ~ "Mixed",
                                 TRUE ~ ClusterType))

sc_vec <- setNames(cluster_type$ClusterType, cluster_type$ID)

d_m_stats$SampleCluster_int <- d_m_stats$SampleCluster

d_m_stats$SampleCluster <- sc_vec[as.character(d_m_stats$ID)]
d_m_stats$SampleCluster <- factor(d_m_stats$SampleCluster, c("AC", "MES/Hypoxia", "GPC", "OPC/NPC", "Neuron", "Mixed"))

d_m_stats %>%
  filter(State %ni% c("Unresolved")) %>%
  ggplot(aes(x = ID, y = State, fill = Freq)) +
  facet_grid(cols = vars(SampleCluster), scales = "free", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_distiller(name = "Proportion", palette = "Purples", direction = 1,
                       breaks = seq(0, .3, .05), limits = c(0, .3), oob = squish, guide = "coloursteps") +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 12, angle = 90))

mal_cluster_stats <- d_m_stats

mal_cluster_data <- mal_cluster_stats %>%
  group_by(Sample, ID, Patient, Timepoint) %>%
  summarise(MalCluster = first(SampleCluster), MalCluster_int = first(SampleCluster_int),
            SI_Mal = shannon_index(Freq),
            .groups = "drop")

gbm_subtypes_tbl <- comp_cluster_data %>%
  left_join(mal_cluster_data %>%
              select(ID, MalCluster, MalCluster_int, SI_Mal), by = "ID")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 5d - Multi-layer graph
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
# Basic pairs (Comp-Mal, Comp-SCP, Mal-SCP)
####################################################################################################################################

feature_scores_all_per_sample$Diff <- apply(feature_scores_all_per_sample[, c("C1", "C3", "C5")] %>% as.matrix(), 1, function(x) sort(x)[3] - sort(x)[2])
feature_scores_all_per_sample$MaxScore <- apply(feature_scores_all_per_sample[, c("C1", "C3", "C5")] %>% as.matrix(), 1, max)

feature_scores_all_per_sample$SCP[feature_scores_all_per_sample$MaxScore < 0 | feature_scores_all_per_sample$Diff < .25] <- "Mixed"

feature_scores_all_per_sample$SCP <- factor(feature_scores_all_per_sample$SCP, c("SCP-ECM", "SCP-Neuronal", "SCP-Glial", "Mixed"))

gbm_subtypes_tbl$SCP <- setNames(feature_scores_all_per_sample$SCP, feature_scores_all_per_sample$Sample)[gbm_subtypes_tbl$Sample]

d <- gbm_subtypes_tbl %>%
  filter(!is.na(MalCluster))

d$Pair1 <- paste0(d$CompCluster, "_", d$MalCluster)
d$Pair2 <- paste0(d$CompCluster, "_", d$SCP)
d$Pair3 <- paste0(d$MalCluster, "_", d$SCP)

pair1_obs_vs_exp <- lapply(unique(d$Pair1), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  nobs <- sum(d$Pair1 == p)
  prob <- (sum(d$CompCluster == p1) / nrow(d)) * (sum(d$MalCluster == p2) / nrow(d))
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / nrow(d), fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = nrow(d), p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair1_obs_vs_exp <- do.call(rbind, pair1_obs_vs_exp)

pair2_obs_vs_exp <- lapply(unique(d$Pair2), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  nobs <- sum(d$Pair2 == p)
  prob <- (sum(d$CompCluster == p1) / nrow(d)) * (sum(d$SCP == p2) / nrow(d))
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / nrow(d), fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = nrow(d), p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair2_obs_vs_exp <- do.call(rbind, pair2_obs_vs_exp)

pair3_obs_vs_exp <- lapply(unique(d$Pair3), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  nobs <- sum(d$Pair3 == p)
  prob <- (sum(d$MalCluster == p1) / nrow(d)) * (sum(d$SCP == p2) / nrow(d))
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / nrow(d), fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = nrow(d), p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair3_obs_vs_exp <- do.call(rbind, pair3_obs_vs_exp)

sig_pairs_tbl <- rbind(pair1_obs_vs_exp %>%
                         filter(pval < .1),
                       pair2_obs_vs_exp %>%
                         filter(pval < .05),
                       pair3_obs_vs_exp %>%
                         filter(pval < .05))

####################################################################################################################################
# Add CNA drivers
####################################################################################################################################

driver_data <- readRDS(paste0(DATA_ROOT, "driver_cna_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -Patient, -Timepoint)

drivers <- colnames(driver_data)[2:6]

driver_data <- driver_data %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, CompCluster, MalCluster, SCP), by = "ID")

d <- driver_data %>%
  filter(!is.na(MalCluster), !is.na(SCP))

N <- length(unique(d$ID))

d <- melt(d, variable.name = "Gene", value.name = "CNA") %>%
  as_tibble() %>%
  filter(CNA == 1) %>%
  mutate(Gene = as.character(Gene))

d$Pair1 <- paste0(d$Gene, "_", d$CompCluster)
d$Pair2 <- paste0(d$Gene, "_", d$MalCluster)
d$Pair3 <- paste0(d$Gene, "_", d$SCP)

pair1_obs_vs_exp <- lapply(unique(d$Pair1), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% driver_data$ID, CompCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair1 == p & d$CNA == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair1_obs_vs_exp <- do.call(rbind, pair1_obs_vs_exp)

pair2_obs_vs_exp <- lapply(unique(d$Pair2), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% driver_data$ID, MalCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair2 == p & d$CNA == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair2_obs_vs_exp <- do.call(rbind, pair2_obs_vs_exp)

pair3_obs_vs_exp <- lapply(unique(d$Pair3), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% driver_data$ID, SCP == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair3 == p & d$CNA == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair3_obs_vs_exp <- do.call(rbind, pair3_obs_vs_exp)

sig_pairs_tbl <- rbind(sig_pairs_tbl,
                       rbind(pair1_obs_vs_exp %>%
                               filter(pval < .1),
                             pair2_obs_vs_exp %>%
                               filter(pval < .1),
                             pair3_obs_vs_exp %>%
                               filter(pval < .1)))

####################################################################################################################################
# Add SNVs
####################################################################################################################################

snv_data <- readRDS(paste0(DATA_ROOT, "driver_snv_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -case_barcode)

snvs <- colnames(snv_data)[2:10]

snv_data <- snv_data %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, CompCluster, MalCluster, SCP), by = "ID")

d <- snv_data %>%
  filter(!is.na(MalCluster), !is.na(SCP))
dim(d)

N <- length(unique(d$ID))

d <- melt(d, variable.name = "Gene", value.name = "SNV") %>%
  as_tibble() %>%
  filter(SNV == 1) %>%
  mutate(Gene = as.character(Gene))

d$Pair1 <- paste0(d$Gene, "_", d$CompCluster)
d$Pair2 <- paste0(d$Gene, "_", d$MalCluster)
d$Pair3 <- paste0(d$Gene, "_", d$SCP)

pair1_obs_vs_exp <- lapply(unique(d$Pair1), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% snv_data$ID, CompCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair1 == p & d$SNV == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair1_obs_vs_exp <- do.call(rbind, pair1_obs_vs_exp)

pair2_obs_vs_exp <- lapply(unique(d$Pair2), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% snv_data$ID, MalCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair2 == p & d$SNV == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair2_obs_vs_exp <- do.call(rbind, pair2_obs_vs_exp)

pair3_obs_vs_exp <- lapply(unique(d$Pair3), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% snv_data$ID, SCP == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair3 == p & d$SNV == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair3_obs_vs_exp <- do.call(rbind, pair3_obs_vs_exp)

sig_pairs_tbl <- rbind(sig_pairs_tbl,
                       rbind(pair1_obs_vs_exp %>%
                               filter(pval < .05),
                             pair2_obs_vs_exp %>%
                               filter(pval < .05),
                             pair3_obs_vs_exp %>%
                               filter(pval < .05)))

####################################################################################################################################
# Generate graph and plot
####################################################################################################################################

edges <- lapply(sig_pairs_tbl$Pair, function(h) {
  
  h <- strsplit(x = as.character(h), split = "_")
  
  vs1 <- h[[1]][1]
  vs2 <- h[[1]][2]
  data.frame(S1 = vs1, S2 = vs2)
})
edges <- do.call(rbind, edges)
edges$log2FC <- sig_pairs_tbl$log2FC
edges$FC <- sig_pairs_tbl$FC

pairs_graph <- graph_from_data_frame(edges, directed = F)

colors <- RColorBrewer::brewer.pal(6, "Spectral")

E(pairs_graph)$weight <- edges$FC
E(pairs_graph)$arrow.size <- edges$FC
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% levels(gbm_subtypes_tbl$SCP)] <- colors[1]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% levels(gbm_subtypes_tbl$MalCluster)] <- colors[2]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% levels(gbm_subtypes_tbl$CompCluster)] <- colors[3]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% drivers] <- colors[4]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% snvs] <- colors[5]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% cnas] <- colors[6]

plot(pairs_graph,
     edge.width = E(pairs_graph)$weight,
     vertex.color = V(pairs_graph)$vertex.color, vertex.label.color="black",
     edge.label.cex = 1,
     edge.curved = F,
     vertex.size = 10, vertex.label.cex = 1.5)

multi_layer_groups <- edges

####################################################################################################################################
# Classify samples to ecosystems
####################################################################################################################################

ecosystems_vec <- components(pairs_graph)$membership
ecosystems_vec <- setNames(paste0("ES", ecosystems_vec), names(ecosystems_vec))

d <- gbm_subtypes_tbl %>%
  select(Sample, ID, CompCluster, MalCluster, SCP) %>%
  mutate(ID = as.character(ID),
         CompCluster = as.character(CompCluster),
         MalCluster = as.character(MalCluster),
         SCP = as.character(SCP))

d <- d %>%
  filter(!is.na(MalCluster), !is.na(SCP))

d <- d %>%
  melt(id.vars = c("Sample", "ID")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable))

d$ES_comp <- ecosystems_vec[d$value]

d <- d %>%
  filter(!is.na(ES_comp))

d <- d %>%
  group_by(Sample, ID, ES_comp) %>%
  mutate(ID = as.character(ID)) %>%
  summarise(n = n(), .groups = "drop")

es_class <- lapply(unique(d$ID), function(id) {
  x <- d %>%
    filter(ID == id) %>%
    arrange(desc(n))
  
  es <- NA
  
  if(nrow(x) == 1 & max(x$n) > 1) {
    es <- x$ES_comp
  } else if(nrow(x) > 1) {
    if(x$n[1] > 2 & x$n[2] <= 1) {
      es <- x$ES_comp[1]
    } else if((x$n[1] == 2 & x$n[2] == 0) | (x$ES_comp[1] == "ES2" & x$n[1] == 2))  {
      es <- x$ES_comp[1]
    }
  }
  
  return(tibble(Sample = x$Sample[1], ID = x$ID[1], ES = es, n = max(x$n)))
})
es_class <- do.call(rbind, es_class)

es_vec <- setNames(rep(NA, nrow(gbm_subtypes_tbl)), gbm_subtypes_tbl$Sample)

es_vec[es_class$Sample] <- es_class$ES

gbm_subtypes_tbl$ES <- es_vec[gbm_subtypes_tbl$Sample]

d <- gbm_subtypes_tbl %>%
  mutate(ES = case_when(CompCluster == "Mixed" | MalCluster == "Mixed" | SCP == "Mixed" ~ NA,
                        TRUE ~ ES))

####################################################################################################################################
# Statistical test for ecosystems
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  select(Sample, ID, CompCluster, MalCluster, SCP) %>%
  mutate(ID = as.character(ID),
         CompCluster = as.character(CompCluster),
         MalCluster = as.character(MalCluster),
         SCP = as.character(SCP))

d <- d %>%
  filter(!is.na(MalCluster), !is.na(SCP))

d <- d %>%
  melt(id.vars = c("Sample", "ID")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable))

perm_test <- lapply(1:25, function(i) {
  print(i)
  
  d_shuff <- d %>%
    group_by(variable) %>%
    mutate(value = sample(x = value, size = n())) %>%
    ungroup()
  
  d_shuff$ES_comp <- ecosystems_vec[d_shuff$value]
  
  d_shuff <- d_shuff %>%
    filter(!is.na(ES_comp))
  
  d_shuff <- d_shuff %>%
    group_by(Sample, ID, ES_comp) %>%
    mutate(ID = as.character(ID)) %>%
    summarise(n = n(), .groups = "drop")
  
  es_class <- lapply(unique(d_shuff$ID), function(id) {
    x <- d_shuff %>%
      filter(ID == id) %>%
      arrange(desc(n))
    
    es <- NA
    
    if(nrow(x) == 1 & max(x$n) > 1) {
      es <- x$ES_comp
    } else if(nrow(x) > 1) {
      if(x$n[1] > 2 & x$n[2] <= 1) {
        es <- x$ES_comp[1]
      } else if((x$n[1] == 2 & x$n[2] == 0) | (x$ES_comp[1] == "ES2" & x$n[1] == 2))  {
        es <- x$ES_comp[1]
      }
    }
    
    return(tibble(Sample = x$Sample[1], ID = x$ID[1], ES = es, n = max(x$n)))
  })
  es_class <- do.call(rbind, es_class)
  
  es_vec <- setNames(rep(NA, nrow(gbm_subtypes_tbl)), gbm_subtypes_tbl$Sample)
  
  es_vec[es_class$Sample] <- es_class$ES
  as_tibble(table(es_vec))
})
perm_test <- do.call(rbind, perm_test)

perm_test_sum <- perm_test %>%
  group_by(es_vec) %>%
  summarise(n = mean(n))

fisher.test(matrix(data = c(perm_test_sum$n, table(gbm_subtypes_tbl$ES)), nrow = 2, byrow = T))

fisher.test(matrix(data = c(c(sum(perm_test_sum$n), 120 - sum(perm_test_sum$n)),
                            c(sum(table(gbm_subtypes_tbl$ES)), 120 - sum(table(gbm_subtypes_tbl$ES)))), nrow = 2, byrow = T))

N <- gbm_subtypes_tbl %>%
  filter(!is.na(MalCluster), !is.na(SCP)) %>%
  nrow()
nobs <- table(gbm_subtypes_tbl$ES)
nexp <- round(perm_test_sum$n, 0)

sapply(1:3, function(i) binom.test(x = nobs[i], n = N, p = nexp[i] / N)$p.value)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 6 - Comprehensive map
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure 6A - ES-classified samples
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

dt_malignant_composition_all_samples <- dT_state(mdata)
dt_overall_composition_all_samples <- dTME_celltype(meta_data)
dt_cell_cycle_all_samples <- dT_CC(mdata)

d <- feature_scores_all_per_sample %>%
  mutate(ID = paste0(Patient, Timepoint))

####################################################################################################################################
# Combine T1 and T2 proportions
####################################################################################################################################

m1 <- cbind(dt_overall_composition_all_samples$T1_ABS, dt_malignant_composition_all_samples$T1_ABS, dt_cell_cycle_all_samples$T1_ABS)
m2 <- cbind(dt_overall_composition_all_samples$T2_ABS, dt_malignant_composition_all_samples$T2_ABS, dt_cell_cycle_all_samples$T2_ABS)

rownames(m1) <- paste0(rownames(m1), "T1")
rownames(m2) <- paste0(rownames(m2), "T2")

m <- rbind(m1, m2)
colnames(m) <- sapply(strsplit(colnames(m), split = "T1_|T2_"), function(x) x[[2]][1])

####################################################################################################################################
# Combine SCPs states
####################################################################################################################################

dm <- melt(d %>%
             select(ID, C1, C3, C5))
dm <- acast(dm, formula = ID ~ variable, value.var = "value")
colnames(dm) <- c("SCP_ECM", "SCP_Neuronal", "SCP_Glial")

dm <- dm[rownames(m), ]

m <- cbind(m, dm)

####################################################################################################################################
# Combine gene-level CNA data
####################################################################################################################################

driver_data <- readRDS(paste0(DATA_ROOT, "driver_cna_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -Patient, -Timepoint)

driver_m <- as.matrix(driver_data[, -1])
rownames(driver_m) <- driver_data$ID

diff_pts <- setdiff(levels(gbm_subtypes_tbl$ID), unique(driver_data$ID))

diff_pts_m <- matrix(data = NA, nrow = length(diff_pts), ncol = ncol(driver_m))
rownames(diff_pts_m) <- diff_pts
colnames(diff_pts_m) <- colnames(driver_m)

driver_m <- rbind(driver_m, diff_pts_m)

####################################################################################################################################
# Combine SNV data
####################################################################################################################################

snv_data <- readRDS(paste0(DATA_ROOT, "driver_snv_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -case_barcode)

snv_m <- as.matrix(snv_data[, -1])
rownames(snv_m) <- snv_data$ID

diff_pts <- setdiff(levels(gbm_subtypes_tbl$ID), unique(snv_data$ID))

diff_pts_m <- matrix(data = NA, nrow = length(diff_pts), ncol = ncol(snv_m))
rownames(diff_pts_m) <- diff_pts
colnames(diff_pts_m) <- colnames(snv_m)

snv_m <- rbind(snv_m, diff_pts_m)

####################################################################################################################################
# Cluster and plot
####################################################################################################################################

d <- d %>%
  left_join(comp_cluster_data %>%
              select(ID, CompCluster, SI_Comp), by = "ID") %>%
  left_join(mal_cluster_data %>%
              select(ID, MalCluster, SI_Mal), by = "ID") %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, ES), by = "ID")

d <- d %>%
  filter(!is.na(ES)) %>%
  filter(Timepoint != "T3")

m_dist <- d %>%
  select(CompCluster, MalCluster, SCP,
         C1, C3, C5) %>%
  cbind(m[d$ID, grep("ABS_TME_", colnames(m))]) %>%
  cbind(m[d$ID, grep("ABS_State_|ABS_CC", colnames(m))])

d_dist <- cluster::daisy(m_dist, metric = "gower")

hc <- fastcluster::hclust(d_dist, method = "average")

d$ID <- factor(d$ID, d$ID[hc$order])

d$ES <- factor(as.character(d$ES), c("ES3", "ES2", "ES1"))

####################################################################################################################################
# Malignant cluster panel
####################################################################################################################################

mal_cluster_color_vec <- setNames(RColorBrewer::brewer.pal(length(levels(gbm_subtypes_tbl$MalCluster)), "Set2"), levels(gbm_subtypes_tbl$MalCluster))

p1 <-
  ggplot(d, aes(x = ID, y = "Malignant", fill = MalCluster)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "MC", values = mal_cluster_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# Composition cluster panel
####################################################################################################################################

comp_cluster_color_vec <- setNames(RColorBrewer::brewer.pal(length(levels(gbm_subtypes_tbl$CompCluster)), "Set1"), levels(gbm_subtypes_tbl$CompCluster))

p2 <-
  ggplot(d, aes(x = ID, y = "Composition", fill = CompCluster)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "CC", values = comp_cluster_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# SCP score panel
####################################################################################################################################

dm <- melt(d %>%
             select(ID, C1, C3, C5, ES),
           id.vars = c("ID", "ES")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable),
         variable = case_when(variable == "C1" ~ "SCP-ECM",
                              variable == "C3" ~ "SCP-Neuronal",
                              variable == "C5" ~ "SCP-Glial"),
         variable = factor(variable, c("SCP-Glial", "SCP-ECM", "SCP-Neuronal")))

p3 <- 
  ggplot(dm, aes(x = ID, y = variable, fill = value)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Score", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-1, 1), oob = squish, labels = c("-1", "", "0", "", "1")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# TME proportion panel
####################################################################################################################################

dm <- melt(m[as.character(d$ID), grep("ABS_TME_", colnames(m))]) %>%
  as_tibble() %>%
  filter(Var2 %ni% grep("_TAMs", as.character(Var2), value = T))

dm <- dm %>%
  group_by(Var1) %>%
  summarise(Malignant = value[Var2 == "ABS_TME_Malignant"], TAM = value[Var2 == "ABS_TME_Macrophage"],
            Oligodendrocyte = value[Var2 == "ABS_TME_Oligodendrocyte"], Astrocyte = value[Var2 == "ABS_TME_Astrocyte"],
            Neuron = value[Var2 == "ABS_TME_Excitatory neuron"] + value[Var2 == "ABS_TME_Inhibitory neuron"],
            Endovascular = value[Var2 == "ABS_TME_Endothel"] + value[Var2 == "ABS_TME_Pericyte"]) %>%
  melt() %>%
  as_tibble() %>%
  mutate(variable = factor(as.character(variable), c("Malignant", "TAM", "Oligodendrocyte", "Astrocyte", "Neuron", "Endovascular")))

dm$Var1 <- factor(dm$Var1, levels(d$ID))

dm <- dm %>%
  group_by(variable) %>%
  mutate(value_c = value - mean(value))

dm <- dm %>%
  left_join(d %>%
              select(ID, ES) %>%
              rename(Var1 = ID), by = "Var1")

p4 <-
  ggplot(dm, aes(x = Var1, y = variable, fill = value_c)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Centered prop", low = "dodgerblue", mid = "white", high = "red",
                       breaks = seq(-.1, .1, .05),
                       limits = c(-.1, .1), oob = squish, labels = c("-.1", "", "0", "", ".1")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# Malignant state proportion panel
####################################################################################################################################

dm <- melt(m[as.character(d$ID), grep("ABS_State_|ABS_CC", colnames(m))]) %>%
  as_tibble() #%>%

dm <- dm %>%
  group_by(Var1) %>%
  summarise(Neuron = value[Var2 == "ABS_State_Neuron"],
            GPC = value[Var2 == "ABS_State_GPC"],
            MES_Hypoxia = value[Var2 == "ABS_State_Hypoxia"] + value[Var2 == "ABS_State_MES"],
            OPC_NPC = value[Var2 == "ABS_State_OPC"] + value[Var2 == "ABS_State_NPC"]) %>%
  melt() %>%
  as_tibble() %>%
  mutate(variable = factor(as.character(variable), c("MES_Hypoxia", "GPC", "OPC_NPC", "Neuron")))

dm$Var1 <- factor(dm$Var1, levels(d$ID))

dm <- dm %>%
  left_join(d %>%
              select(ID, ES) %>%
              rename(Var1 = ID), by = "Var1")

p5 <-
  ggplot(dm, aes(x = Var1, y = variable, fill = value)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_distiller(name = "Prop", palette = "Purples", direction = 1,
                       breaks = seq(0, .3, .05), limits = c(0, .3), oob = squish, guide = "coloursteps") +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# Driver CNA panel
####################################################################################################################################
dm <- melt(driver_m) %>%
  as_tibble() %>%
  filter(Var2 %in% c("EGFR.amp", "MDM2.amp", "CDKN2A.del"))

dm <- dm %>%
  filter(as.character(Var1) %in% levels(d$ID))

dm$Var1 <- factor(dm$Var1, levels(d$ID))

dm$value <- as.character(dm$value)
table(dm$value)

dm$value[is.na(dm$value)] <- "NA"

dm <- dm %>%
  left_join(d %>%
              select(ID, ES) %>%
              rename(Var1 = ID), by = "Var1")

p6 <-
  ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "", values = c("1" = "black", "0" = "white", "NA" = "grey"),
                    labels = c("0" = "Not detected", "1" = "Detected", "NA" = "NA")) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# SNV panel
####################################################################################################################################

dm <- melt(snv_m[, c("TP53", "NF1", "RB1")]) %>%
  as_tibble()

dm <- dm %>%
  filter(as.character(Var1) %in% levels(d$ID))

dm$Var1 <- factor(as.character(dm$Var1), levels(d$ID))
dm$Var2 <- factor(dm$Var2, c("NF1", "TP53", "RB1"))

dm$value <- as.character(dm$value)
table(dm$value)

dm$value[is.na(dm$value)] <- "NA"

dm <- dm %>%
  left_join(d %>%
              select(ID, ES) %>%
              rename(Var1 = ID), by = "Var1")

p7 <-
  ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  facet_grid(cols = vars(ES), scales = "free_x", space = "free_x") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "", values = c("1" = "black", "0" = "white", "NA" = "grey"),
                    labels = c("0" = "Not detected", "1" = "Detected", "NA" = "NA")) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_text(size = 10, angle = 90),
                 axis.text.y = element_text(size = 12),
                 legend.position = "right")

p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 1, heights = c(.5, .5, 2, 3, 2, 1, 1))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Figure EDF10d - ES-unclassified samples
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

d <- feature_scores_all_per_sample %>%
  mutate(ID = paste0(Patient, Timepoint))

####################################################################################################################################
# Combine T1 and T2 proportions
####################################################################################################################################

m1 <- cbind(dt_overall_composition_all_samples$T1_ABS, dt_malignant_composition_all_samples$T1_ABS, dt_cell_cycle_all_samples$T1_ABS)
m2 <- cbind(dt_overall_composition_all_samples$T2_ABS, dt_malignant_composition_all_samples$T2_ABS, dt_cell_cycle_all_samples$T2_ABS)

rownames(m1) <- paste0(rownames(m1), "T1")
rownames(m2) <- paste0(rownames(m2), "T2")

m <- rbind(m1, m2)
colnames(m) <- sapply(strsplit(colnames(m), split = "T1_|T2_"), function(x) x[[2]][1])

####################################################################################################################################
# Combine SCPs states
####################################################################################################################################

dm <- melt(d %>%
             select(ID, C1, C3, C5))
dm <- acast(dm, formula = ID ~ variable, value.var = "value")
colnames(dm) <- c("SCP_ECM", "SCP_Neuronal", "SCP_Glial")

dm <- dm[rownames(m), ]

m <- cbind(m, dm)

####################################################################################################################################
# Combine gene-level CNA data
####################################################################################################################################

driver_data <- readRDS(paste0(DATA_ROOT, "driver_cna_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -Patient, -Timepoint)

driver_m <- as.matrix(driver_data[, -1])
rownames(driver_m) <- driver_data$ID

diff_pts <- setdiff(levels(gbm_subtypes_tbl$ID), unique(driver_data$ID))

diff_pts_m <- matrix(data = NA, nrow = length(diff_pts), ncol = ncol(driver_m))
rownames(diff_pts_m) <- diff_pts
colnames(diff_pts_m) <- colnames(driver_m)

driver_m <- rbind(driver_m, diff_pts_m)

####################################################################################################################################
# Combine SNV data
####################################################################################################################################

snv_data <- readRDS(paste0(DATA_ROOT, "driver_snv_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -case_barcode)

snv_m <- as.matrix(snv_data[, -1])
rownames(snv_m) <- snv_data$ID

diff_pts <- setdiff(levels(gbm_subtypes_tbl$ID), unique(snv_data$ID))

diff_pts_m <- matrix(data = NA, nrow = length(diff_pts), ncol = ncol(snv_m))
rownames(diff_pts_m) <- diff_pts
colnames(diff_pts_m) <- colnames(snv_m)

snv_m <- rbind(snv_m, diff_pts_m)

####################################################################################################################################
# Cluster and plot
####################################################################################################################################

d <- d %>%
  left_join(comp_cluster_data %>%
              select(ID, CompCluster, SI_Comp), by = "ID") %>%
  left_join(mal_cluster_data %>%
              select(ID, MalCluster, SI_Mal), by = "ID") %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, ES), by = "ID")

d <- d %>%
  filter(is.na(ES)) %>%
  filter(Timepoint != "T3")

m_dist <- d %>%
  select(CompCluster, MalCluster, SCP,
         C1, C3, C5) %>%
  cbind(m[d$ID, grep("ABS_TME_", colnames(m))]) %>%
  cbind(m[d$ID, grep("ABS_State_|ABS_CC", colnames(m))])

d_dist <- cluster::daisy(m_dist, metric = "gower")

hc <- fastcluster::hclust(d_dist, method = "average")

d$ID <- factor(d$ID, d$ID[hc$order])

####################################################################################################################################
# Malignant cluster panel
####################################################################################################################################

mal_cluster_color_vec <- setNames(RColorBrewer::brewer.pal(length(levels(gbm_subtypes_tbl$MalCluster)), "Set2"), levels(gbm_subtypes_tbl$MalCluster))

p1 <-
  ggplot(d, aes(x = ID, y = "Malignant", fill = MalCluster)) +
  geom_tile(color = "black") +
  scale_fill_manual(name = "MC", values = mal_cluster_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# Composition cluster panel
####################################################################################################################################

comp_cluster_color_vec <- setNames(RColorBrewer::brewer.pal(length(levels(gbm_subtypes_tbl$CompCluster)), "Set1"), levels(gbm_subtypes_tbl$CompCluster))

p2 <-
  ggplot(d, aes(x = ID, y = "Composition", fill = CompCluster)) +
  geom_tile(color = "black") +
  scale_fill_manual(name = "CC", values = comp_cluster_color_vec) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# SCP score panel
####################################################################################################################################

dm <- melt(d %>%
             select(ID, C1, C3, C5),
           id.vars = c("ID")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable),
         variable = case_when(variable == "C1" ~ "SCP-ECM",
                              variable == "C3" ~ "SCP-Neuronal",
                              variable == "C5" ~ "SCP-Glial"),
         variable = factor(variable, c("SCP-ECM", "SCP-Glial", "SCP-Neuronal")))

p3 <- 
  ggplot(dm, aes(x = ID, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Score", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-1, 1), oob = squish, labels = c("-1", "", "0", "", "1")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# TME proportion panel
####################################################################################################################################

dm <- melt(m[as.character(d$ID), grep("ABS_TME_", colnames(m))]) %>%
  as_tibble() %>%
  filter(Var2 %ni% grep("_TAMs", as.character(Var2), value = T))

dm <- dm %>%
  group_by(Var1) %>%
  summarise(Malignant = value[Var2 == "ABS_TME_Malignant"], TAM = value[Var2 == "ABS_TME_Macrophage"],
            Oligodendrocyte = value[Var2 == "ABS_TME_Oligodendrocyte"], Astrocyte = value[Var2 == "ABS_TME_Astrocyte"],
            Neuron = value[Var2 == "ABS_TME_Excitatory neuron"] + value[Var2 == "ABS_TME_Inhibitory neuron"],
            Endovascular = value[Var2 == "ABS_TME_Endothel"] + value[Var2 == "ABS_TME_Pericyte"]) %>%
  melt() %>%
  as_tibble() %>%
  mutate(variable = factor(as.character(variable), c("Malignant", "TAM", "Oligodendrocyte", "Astrocyte", "Neuron", "Endovascular")))

dm$Var1 <- factor(dm$Var1, levels(d$ID))

dm <- dm %>%
  group_by(variable) %>%
  mutate(value_c = value - mean(value))

p4 <-
  ggplot(dm, aes(x = Var1, y = variable, fill = value_c)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Centered prop", low = "dodgerblue", mid = "white", high = "red",
                       breaks = seq(-.1, .1, .05),
                       limits = c(-.1, .1), oob = squish, labels = c("-.1", "", "0", "", ".1")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# Malignant state proportion panel
####################################################################################################################################

dm <- melt(m[as.character(d$ID), grep("ABS_State_|ABS_CC", colnames(m))]) %>%
  as_tibble() #%>%

dm <- dm %>%
  group_by(Var1) %>%
  summarise(GPC = value[Var2 == "ABS_State_GPC"],
            MES_Hypoxia = value[Var2 == "ABS_State_Hypoxia"] + value[Var2 == "ABS_State_MES"],
            OPC_NPC = value[Var2 == "ABS_State_OPC"] + value[Var2 == "ABS_State_NPC"],
            NEU = value[Var2 == "ABS_State_Neuron"],
            AC_CILIA = value[Var2 == "ABS_State_AC"] + value[Var2 == "ABS_State_Cilia"]) %>%
  melt() %>%
  as_tibble() %>%
  mutate(variable = factor(as.character(variable), c("MES_Hypoxia", "GPC", "AC_CILIA", "OPC_NPC", "NEU")))

dm$Var1 <- factor(dm$Var1, levels(d$ID))

p5 <-
  ggplot(dm, aes(x = Var1, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_distiller(name = "Prop", palette = "Purples", direction = 1,
                       breaks = seq(0, .3, .05), limits = c(0, .3), oob = squish, guide = "coloursteps") +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# Driver CNA panel
####################################################################################################################################

dm <- melt(driver_m) %>%
  as_tibble() %>%
  filter(Var2 %in% c("EGFR.amp", "MDM2.amp", "CDKN2A.del"))

dm <- dm %>%
  filter(as.character(Var1) %in% levels(d$ID))

dm$Var1 <- factor(dm$Var1, levels(d$ID))

dm$value <- as.character(dm$value)

dm$value[is.na(dm$value)] <- "NA"

p6 <-
  ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_manual(name = "", values = c("1" = "black", "0" = "white", "NA" = "grey"),
                    labels = c("0" = "Not detected", "1" = "Detected", "NA" = "NA")) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12),
                 legend.position = "right")

####################################################################################################################################
# SNV panel
####################################################################################################################################

dm <- melt(snv_m[, c("TP53", "NF1", "RB1")]) %>%
  as_tibble()

dm <- dm %>%
  filter(as.character(Var1) %in% levels(d$ID))

dm$Var1 <- factor(as.character(dm$Var1), levels(d$ID))
dm$Var2 <- factor(dm$Var2, c("NF1", "TP53", "RB1"))

dm$value <- as.character(dm$value)
table(dm$value)

dm$value[is.na(dm$value)] <- "NA"

p7 <-
  ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_manual(name = "", values = c("1" = "black", "0" = "white", "NA" = "grey"),
                    labels = c("0" = "Not detected", "1" = "Detected", "NA" = "NA")) +
  xlab("") +
  ylab("") +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_text(size = 10, angle = 90),
                 axis.text.y = element_text(size = 12),
                 legend.position = "right")

p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 1, heights = c(.5, .5, 1, 3, 2, 1, 1))
