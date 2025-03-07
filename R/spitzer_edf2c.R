#################################
## Title: Extended Data Figure 2c in Spitzer et al - summary of composition, malignant states and baseline profile groups
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce heatmaps summarizing the composition, malignant states and baseline profile groups
#################################

####################################################################################################################################
# Composition panel
####################################################################################################################################

comp_cluster_stats <- comp_cluster_stats %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, PvsR), by = "ID")

comp_cluster_stats_summarized <- comp_cluster_stats %>%
  group_by(SampleCluster, CellType) %>%
  summarise(Mean = mean(Freq), .groups = "drop") %>%
  mutate(CellType = factor(as.character(CellType), c("Malignant", "TAM", "Oligodendrocyte", "Astrocyte", "Excitatory neuron", "Inhibitory neuron",
                                                     "Endothel", "Pericyte", "Lymphocyte", "OPC", "Other")))

comp_cluster_stats_summarized <- melt(comp_cluster_stats_summarized) %>%
  as_tibble()

comp_cluster_stats_summarized %>%
  ggplot(aes(x = SampleCluster, y = CellType, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(name = "Proportion", colours = c("white","#FFFFCC", "#FFEDA0", "#FED976",
                                                        "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C",
                                                        "#BD0026", "#800026"), limits = c(0, .8), oob = squish) +
  xlab("Assigned cluster") +
  ylab("Cell type") +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1, ticks.colour = "black")) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 16))

####################################################################################################################################
# Malignant state panel
####################################################################################################################################

mal_cluster_stats <- mal_cluster_stats %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, PvsR), by = "ID")

mal_cluster_stats_summarized <- mal_cluster_stats %>%
  group_by(SampleCluster, State) %>%
  summarise(Mean = mean(Freq), .groups = "drop")

mal_cluster_stats_summarized <- melt(mal_cluster_stats_summarized) %>%
  as_tibble()

mal_cluster_stats_summarized %>%
  ggplot(aes(x = SampleCluster, y = State, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_distiller(name = "Proportion", palette = "Purples", direction = 1,
                       breaks = seq(0, .3, .05), limits = c(0, .3), oob = squish, guide = "coloursteps") +
  xlab("Assigned cluster") +
  ylab("Malignant state") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 16))

####################################################################################################################################
# Baseline profile panel
####################################################################################################################################

gbm_subtypes_tbl <- gbm_subtypes_tbl %>%
  left_join(sample_data %>%
              select(Sample, PvsR), by = "Sample") %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrence"))
table(gbm_subtypes_tbl$PvsR)

dm <- melt(feature_scores_all_per_sample %>%
             mutate(ID = paste0(Patient, Timepoint)) %>%
             select(ID, C1, C3, C5, PvsR, SCP),
           id.vars = c("ID", "PvsR", "SCP")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable),
         variable = case_when(variable == "C1" ~ "SCP-ECM",
                              variable == "C3" ~ "SCP-Neuronal",
                              variable == "C5" ~ "SCP-Glial"),
         variable = factor(variable, c("SCP-ECM", "SCP-Neuronal", "SCP-Glial")))

dm_m <- acast(dm, formula = ID ~ variable, value.var = "value")

hc <- fastcluster::hclust(d = as.dist(1 - cor(t(dm_m))), method = "average")

dm$ID <- factor(dm$ID, hc$labels[hc$order])

scp_stats_summarized <- dm %>%
  group_by(variable, SCP) %>%
  summarise(Mean = mean(value), .groups = "drop")

scp_stats_summarized %>%
  ggplot(aes(x = SCP, y = variable, fill = Mean)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Score", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-1, 1), oob = squish,
                       breaks = c(-1, -.5, 0, .5, 1), labels = c("-1", "", "0", "", "1")) +
  xlab("Assigned SCP") +
  ylab("Baseline program") +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), axis.text.x = element_text(size = 16))
