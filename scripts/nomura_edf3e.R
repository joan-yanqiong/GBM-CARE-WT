#################################
## Title: Extended Data Figure 3 panel e in Nomura et al - similarity with states from previous GBM studies
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a heatmap showing similarity between CARE states with GBM studies other than Neftel et al (Bhaduri et al, Mathur et al)
#################################

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
