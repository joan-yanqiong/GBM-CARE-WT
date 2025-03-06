#################################
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