#################################
## Title: Extended Data Figure 2 panel f in Nomura et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Pie chart demonstrating the proportion (%) of each cell type in the Smart-seq2 dataset
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

load("data/CARE_SS2_allProtocols.RData", verbose = T)
mData <- aSeurat@meta.data %>%
  as_tibble
umap <- as.data.frame(aSeurat@reductions$umap@cell.embeddings)

ccol <- c("Endothelial" = "#FCCDE5", "Lymphocyte" = "#8DD3C7", "Malignant" = "#FB8072", "Neuron" = "#BC80BD", "Oligodendrocyte" = "#B3DE69", "Pericyte" = "#FFFFB3", "TAM" = "#80B1D3", "Unresolved" = "grey")

df <- table(mData$CellType)
df <- data.frame(group = names(df), value = as.numeric(df))
df <- cbind(df, perc = round((df$value/sum(df$value))*100, 1))
df <- cbind(df, labels = paste0(df$group, "\n(", df$perc, "%)"))
df <- df[order(df$perc), ]
df$group <- factor(df$group, levels = df$group)
df <- cbind(df, labelsIn = df$labels, labelsOut = df$labels)
df[1:5, "labelsIn"] <- ""
df[6:8, "labelsOut"] <- ""

df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = ccol[levels(df$group)]) +
  geom_text(aes(label = labelsIn), position = position_stack(vjust = 0.6)) +
  geom_label_repel(data = df2, aes(y = pos, label = labelsOut), nudge_x = 1, show.legend = FALSE, fill = NA, label.size = NA) +
  theme_void() +
  NoLegend()


