#################################
<<<<<<< HEAD
## Title: Extended Data Figure 1a-b in Spitzer et al - number of cells and average number of detected genes across timepoints
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce boxplots of number of cells and average number of detected genes across timepoints
#################################

####################################################################################################################################
# Number of cells per sample across timepoints
####################################################################################################################################

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Ncells = n()) %>%
  mutate(Timepoint = factor(Timepoint, c("T1", "T2", "T3"))) %>%
  ggboxplot(x = "Timepoint", y = "Ncells", add = "jitter", fill = "Timepoint") +
  scale_color_brewer(name = "", palette = "Set1") +
  xlab("Timepoint") +
  ylab("Number of cells\nper sample") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "none") +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Ncells = n()) %>%
  ungroup() %>%
  summarise(p1 = wilcox.test(x = Ncells[Timepoint == "T1"], y = Ncells[Timepoint == "T2"])$p.value,
            p2 = wilcox.test(x = Ncells[Timepoint == "T1"], y = Ncells[Timepoint == "T3"])$p.value,
            p3 = wilcox.test(x = Ncells[Timepoint == "T2"], y = Ncells[Timepoint == "T3"])$p.value)

####################################################################################################################################
# Average number of detected genes per sample across timepoints
####################################################################################################################################

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Comp = mean(nFeature_RNA)) %>%
  mutate(Timepoint = factor(Timepoint, c("T1", "T2", "T3"))) %>%
  ggboxplot(x = "Timepoint", y = "Comp", add = "jitter", fill = "Timepoint") +
  scale_color_brewer(name = "", palette = "Set1") +
  xlab("Timepoint") +
  ylab("Average number of detected\ngenes per sample") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "none") +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Comp = mean(nFeature_RNA)) %>%
  ungroup() %>%
  summarise(p1 = wilcox.test(x = Comp[Timepoint == "T1"], y = Comp[Timepoint == "T2"])$p.value,
            p2 = wilcox.test(x = Comp[Timepoint == "T1"], y = Comp[Timepoint == "T3"])$p.value,
            p3 = wilcox.test(x = Comp[Timepoint == "T2"], y = Comp[Timepoint == "T3"])$p.value)
=======
## Title: Extended Data Figure 1 panel a_b in Spitzer et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Number of cells per sample (post QC) and Average number of detected genes per sample (post QC) across the 3 time-points
#################################

# Necessary packages
library(tidyverse)
library(Matrix)
library(corrplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(colorspace)
library(circlize)
library(Seurat)
library(ggplot2)
library(ggrepel)

load("data/CARE_10x_pData_240722.RData", verbose = T)
mData <- pData

toAdd <- mData %>%
  group_by(ID) %>%
  summarise(nCells = n(),
            Complexity = median(Complexity),
            mitoPerc = median(mitoPerc)
  )

pData <- mData %>%
  distinct(ID, .keep_all = T) %>%
  select(ID, Patient, Timepoint, CompCluster, MalCluster, SCP, ES, TumorLocation, TumorLocationGroup, SiteRecurrence, Radiation, AlkylateAgents, Steroid) %>%
  left_join(toAdd, by = "ID")
rm(toAdd)

colTimepoint <- c("#F3766E", "#2AB34B", "#7094CD")
names(colTimepoint) <- c("T1", "T2", "T3")

pData %>%
  ggplot(aes(x = Timepoint, y = nCells, fill = Timepoint, col = Timepoint)) +
  theme_bw(base_size = 14) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = colTimepoint) +
  scale_fill_manual(values = rep("white", 4)) +
  new_scale_fill() +
  geom_point(aes(fill = Timepoint), pch = 21, size = 1.5, col = "black", position = position_jitterdodge()) +
  scale_fill_manual(values = colTimepoint) +
  labs(title = "", x = "Timepoint", y = "Number of cells per sample") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

pData %>%
  ggplot(aes(x = Timepoint, y = Complexity, fill = Timepoint, col = Timepoint)) +
  theme_bw(base_size = 14) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = colTimepoint) +
  scale_fill_manual(values = rep("white", 4)) +
  new_scale_fill() +
  geom_point(aes(fill = Timepoint), pch = 21, size = 1.5, col = "black", position = position_jitterdodge()) +
  scale_fill_manual(values = colTimepoint) +
  labs(title = "", x = "Timepoint", y = "Average number of detected genes per sample") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
>>>>>>> 0f47ea65bff05825c5ef4c69fc28c985294abdd9
