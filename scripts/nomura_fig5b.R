#################################
## Title: Figure 5 panel b in Nomura et al - heatmap of composition groups
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Generate a heatmap showing the composition groups across all samples
#################################

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
