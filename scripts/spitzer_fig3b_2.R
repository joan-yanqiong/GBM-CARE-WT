#################################
## Title: Spitzer Figure 3 panel b - examples of longitudinal heterogenity
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate heatmaps of features showing examples of longitudinal heterogeneity
#################################

d <- gbm_subtypes_tbl %>%
  group_by(Patient) %>%
  summarise(C1 = CompCluster[Timepoint == "T1"], C2 = CompCluster[Timepoint == "T2"],
            M1 = MalCluster[Timepoint == "T1"], M2 = MalCluster[Timepoint == "T2"],
            CTraj = paste0(C1, "-", C2), MTraj = paste0(M1, "-", M2),
            CTrajOpp = paste0(C2, "-", C1), MTrajOpp = paste0(M2, "-", M1),
            .groups = "drop")

for(i in 1:length(d)) {
  for(j in 1:length(d)) {
    if(i != j & (d$MTraj[i] == d$MTrajOpp[j]))
      print(paste0(i, ", ", j))
  }
}

smpls <- mdata$Sample[mdata$Patient %in% d$Patient[1:2]] %>% unique()

pts <- gbm_subtypes_tbl %>%
  group_by(Patient) %>%
  filter("T1" %in% Timepoint & "T2" %in% Timepoint & "T3" %in% Timepoint) %>%
  filter(Patient != "P58") %>%
  pull(Patient) %>%
  unique()

smpls <- c(smpls, mdata$Sample[mdata$Patient %in% pts] %>% unique())

####################################################################################################################################
# Composition panel
####################################################################################################################################

d_stats <- meta_data %>%
  filter(Sample %in% smpls) %>%
  group_by(Patient, Timepoint, ID_factor, Sample, CellType_reduced) %>%
  summarise(n = n()) %>%
  group_by(Patient, Timepoint, ID_factor, Sample) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup()
d_stats <- acast(d_stats, formula = ID_factor ~ CellType_reduced, value.var = "Freq")
d_stats[is.na(d_stats)] <- 0
d_stats <- melt(d_stats) %>%
  as_tibble()
colnames(d_stats) <- c("ID_factor", "CellType", "Freq")
d_stats <- d_stats %>%
  left_join(sample_data %>%
              select(ID_factor, Patient, Timepoint, Sample), by = "ID_factor")

d_stats <- d_stats %>%
  filter(CellType %in% c("Malignant", "Oligodendrocyte", "TAM", "Astrocyte", "Excitatory neuron")) %>%
  mutate(CellType = factor(CellType, c("Malignant", "TAM", "Oligodendrocyte",
                                       "Astrocyte", "Excitatory neuron")))

d_stats <- d_stats %>%
  mutate(Timepoint = factor(Timepoint, c("T3", "T2", "T1")))

p1 <- 
  ggplot(d_stats, aes(x = CellType, y = Timepoint, fill = Freq)) +
  facet_grid(rows = vars(Patient), scales = "free_y", space = "free_y") +
  geom_tile(color = "black") +
  scale_fill_gradient(name = "Proportion", low = "white", high = "#006D2C",
                      oob = squish, limits = c(0, .8), breaks = seq(0, .8, .1)) +
  xlab("Cell type") +
  ylab("Timepoint") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 90))

####################################################################################################################################
# State panel
####################################################################################################################################

d_stats <- mdata %>%
  filter(Sample %in% smpls) %>%
  group_by(Patient, Timepoint, ID_factor, Sample, State) %>%
  summarise(n = n()) %>%
  group_by(Patient, Timepoint, ID_factor, Sample) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup()
d_stats <- acast(d_stats, formula = ID_factor ~ State, value.var = "Freq")
d_stats[is.na(d_stats)] <- 0
d_stats <- melt(d_stats) %>%
  as_tibble()
colnames(d_stats) <- c("ID_factor", "State", "Freq")
d_stats <- d_stats %>%
  left_join(sample_data %>%
              select(ID_factor, Patient, Timepoint, Sample), by = "ID_factor")

d_stats <- d_stats %>%
  filter(State != "Unresolved")

cc_stats <- mdata %>%
  filter(Sample %in% smpls) %>%
  group_by(Patient, Timepoint, ID_factor, Sample) %>%
  summarise(n = sum(isCC), N = n(), Freq = n / N, .groups = "drop") %>%
  mutate(State = "Cycling") %>%
  select(ID_factor, State, Freq, Patient, Timepoint, Sample)

d_stats <- rbind(d_stats, cc_stats)

d_stats$State <- factor(d_stats$State, c("Cilia", "AC", "MES", "Hypoxia", "Stress",
                                         "GPC", "OPC", "NPC", "Neuron", "Cycling"))

d_stats <- d_stats %>%
  mutate(Timepoint = factor(Timepoint, c("T3", "T2", "T1")))

p2 <- 
  ggplot(d_stats, aes(x = State, y = Timepoint, fill = Freq)) +
  facet_grid(rows = vars(Patient), scales = "free_y", space = "free_y") +
  geom_tile(color = "black") +
  scale_fill_distiller(name = "Proportion", palette = "Purples", direction = 1,
                       breaks = seq(0, .5, .1), limits = c(0, .5), oob = squish, guide = "coloursteps") +
  xlab("State") +
  ylab("") +
  guides(fill = guide_colorsteps(frame.colour = "black", frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 90),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

####################################################################################################################################
# BP panel
####################################################################################################################################

d_stats <- feature_scores_all_per_sample %>%
  filter(Sample %in% smpls) %>%
  select(Patient, Timepoint, Sample, C1, C3, C5) %>%
  melt() %>%
  as_tibble() %>%
  rename(BP = variable, Score = value) %>%
  mutate(BP = as.character(BP)) %>%
  mutate(BP = case_when(BP == "C1" ~ "ECM",
                        BP == "C3" ~ "Neuronal",
                        BP == "C5" ~ "Glial")) %>%
  mutate(BP = factor(BP, c("ECM", "Neuronal", "Glial")))

d_stats <- d_stats %>%
  mutate(Timepoint = factor(Timepoint, c("T3", "T2", "T1")))

p3 <- 
  ggplot(d_stats, aes(x = BP, y = Timepoint, fill = Score)) +
  facet_grid(rows = vars(Patient), scales = "free_y", space = "free_y") +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", midpoint = 0,
                       oob = squish, limits = c(-1.5, 1.5)) +
  xlab("Baseline Profile") +
  ylab("") +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1, ticks.colour = "black")) +
  theme_gbm_pvsr(panel.spacing = unit(.5, "lines"), legend.position = "top") +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 90),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

p1 + p2 + p3 + plot_layout(nrow = 1, widths = c(1, 2, 0.5))
