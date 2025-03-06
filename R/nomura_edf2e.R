#################################
## Title: Extended Data Figure 2 panel e in Nomura et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Heatmaps representing the average gene expression level (log2 relative expression) of cell type marker genes within each TME cell type for the Smart-seq2 dataset.
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

load("data/CARE_SS2_allProtocols.RData", verbose = T)
mData <- aSeurat@meta.data %>%
  as_tibble
umap <- as.data.frame(aSeurat@reductions$umap@cell.embeddings)

aSeurat <- subset(aSeurat, subset = CellType != "Malignant")
aSeurat <- subset(aSeurat, subset = CellType != "Unresolved")
mData <- aSeurat@meta.data
mData$CellType <- factor(mData$CellType, levels = c("Neuron", "Lymphocyte", "Pericyte", "Endothelial", "TAM", "Oligodendrocyte"))
aSeurat@meta.data <- mData

aSeurat1 <- subset(aSeurat, subset = Protocol == "Protocol 1")
aSeurat2 <- subset(aSeurat, subset = Protocol == "Protocol 2")
mData1 <- aSeurat1@meta.data
mData2 <- aSeurat2@meta.data

ffeat <- c("MAG", "MBP", #oligo
           "CD163", "CD74", #mye
           "ENG", "VWF", #endo
           "DCN", "PDGFRB", #peri
           "CD8A", "GZMA", #Tcel
           "SYT1", "SNAP25" #Neuro
)
geData1 <- aSeurat1@assays$RNA@data[ffeat, ]
geData1 <- t(apply(geData1, 1, \(x) tapply(x, mData1$CellType, mean)))
colnames(geData1) <- str_c("P1_", colnames(geData1))
geData2 <- aSeurat2@assays$RNA@data[ffeat, ]
geData2 <- t(apply(geData2, 1, \(x) tapply(x, mData2$CellType, mean)))
colnames(geData2) <- str_c("P2_", colnames(geData2))
geData <- cbind(geData1, geData2)
geData <- geData[, c("P1_Oligodendrocyte", "P1_TAM", "P1_Endothelial", "P1_Pericyte", "P1_Lymphocyte", "P1_Neuron", "P2_Oligodendrocyte", "P2_TAM", "P2_Lymphocyte")]
geData <- t(apply(geData, 1, \(x) (x - mean(x))/sd(x)))
colnames(geData) <- gsub("P1_|P2_", "", colnames(geData))
range(geData)


library(corrplot)
library(RColorBrewer)
par(mar=c(8, 5, 8, 2.1))
corrplot(geData, method="color", is.corr=F,
         pch.cex = 3, pch.col = "white",
         col=colorRampPalette(c("blue4", "white", "firebrick"))(75),
         tl.srt = 45, tl.col = "black", tl.cex = 1, cl.pos = "b", mar = c(0, 0, 0, 5), col.lim = c(-3.1, 3.1))

