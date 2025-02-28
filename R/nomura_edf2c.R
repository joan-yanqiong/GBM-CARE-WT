#################################
## Title: Extended Data Figure 2 panel c in Nomura et al - Heatmap of cell type markers
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Heatmap of cell type markers
#################################

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
