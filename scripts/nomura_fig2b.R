#################################
## Title: Figure 2 panel b in Nomura et al - heatmap of state scores
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a heatmap of per-cell scores for each of the new states
#################################

d <- mdata

# Downsample to 20K cells for easier plotting
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
