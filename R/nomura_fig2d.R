#################################
## Title: Figure 2 panel d in Nomura et al - State similarity with states from normal brain development
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a heatmap showing the similarity between the malignant states and normal brain development states
#################################

LIU_markers <- read.csv(paste0(DATA_ROOT, "liu_marker_genes.csv"), header = T, stringsAsFactors = F)

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
