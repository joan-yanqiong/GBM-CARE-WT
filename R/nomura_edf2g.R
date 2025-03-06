#################################
## Title: Extended Data Figure 2 panel g in Nomura et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Boxplots of cell type proportions (%) in each tumor (n=44) in both 10x (red) and Smart-seq2.
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
mData <- aSeurat@meta.data
mDataSS2 <- mData

IDs <- sort(unique(mData$ID))

pData10x <- readRDS("data/celltype_meta_data_2025_01_08.RDS")
mData <- pData %>%
  filter(ID %in% IDs)

mData[mData$CellType == "Endothel", "CellType"] <- "Endothelial"
mData[mData$CellType %in% c("Excitatory neuron", "Inhibitory neuron"), "CellType"] <- "Neuron"
mData[!mData$CellType %in% c("Malignant", "TAM", "Oligodendrocyte", "Neuron", "Endothelial", "Pericyte", "Lymphocyte", "ParsTuber", "RG"), "CellType"] <- "Other"

mData10x <- mData

mDataSS2 <- table(mDataSS2$ID, mDataSS2$CellType)
mDataSS2 <- as.matrix(mDataSS2)
mDataSS2 <- t(apply(mDataSS2, 1, function(x) (x/sum(x))*100))
mDataSS2 <- mDataSS2[, c("Malignant", "TAM", "Oligodendrocyte", "Neuron", "Endothelial", "Pericyte", "Lymphocyte")]
mDataSS2 <- as.data.frame(mDataSS2)

mData10x <- table(mData10x$ID, mData10x$CellType)
mData10x <- as.matrix(mData10x)
mData10x <- t(apply(mData10x, 1, function(x) (x/sum(x))*100))
mData10x <- mData10x[, c("Malignant", "TAM", "Oligodendrocyte", "Neuron", "Endothelial", "Pericyte", "Lymphocyte")]
mData10x <- as.data.frame(mData10x)

pData <- NULL
for(i in 1:ncol(mData10x)){
  tmp <- tibble(Type = "10x", ID = rownames(mData10x), Patient = ID %>%
                  str_split("T") %>%
                  map_chr(function(x) x[1]), Timepoint = ID %>% str_split("T") %>%
                  map_chr(function(x) x[2]) %>% str_c("T", .),
                CellType = colnames(mData10x)[i], perc = mData10x[, i]
                ) %>%
    bind_rows(tibble(Type = "SS2", ID = rownames(mDataSS2), Patient = ID %>%
                       str_split("T") %>%
                       map_chr(function(x) x[1]),
                     Timepoint = ID %>%
                       str_split("T") %>%
                       map_chr(function(x) x[2]) %>%
                       str_c("T", .),
                     CellType = colnames(mDataSS2)[i], perc = mDataSS2[, i])
              )
  
  pData <- bind_rows(pData, tmp)
}

pData$CellType <- factor(pData$CellType, levels = c("Malignant", "TAM", "Oligodendrocyte", "Neuron", "Endothelial", "Pericyte", "Lymphocyte"))

pData %>%
  ggplot(aes(x = Type, y = perc)) +
  theme_bw(base_size = 14) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Type), color = "black", pch = 21, size = 3) +
  geom_line(aes(group = ID)) +
  facet_wrap(~ CellType, scales = "free_y", nrow = 1) +
  labs(title = "", y = "Proportions") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")
  )
