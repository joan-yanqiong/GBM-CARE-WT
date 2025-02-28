#################################
## Title: Figure 6 panel a in Nomura et al - layers heatmap of samples classified to ecosystems
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Generate a layers heatmap of samples classified to ecosystems
#################################

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
