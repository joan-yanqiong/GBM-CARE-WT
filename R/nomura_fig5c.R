#################################
## Title: Figure 5 panel c in Nomura et al - heatmap of malignant state groups
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Generate a heatmap showing the malignant state groups across all samples
#################################

# Select for the analysis samples with at least 50 malignant cells
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
