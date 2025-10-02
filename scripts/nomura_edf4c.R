#################################
## Title: Extended Data Figure 4 panel c in Nomura et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Pearson’s correlation of the proportions (%) of cells classified in each MP per tumor in both 10x (rows) and Smart-seq2 (columns)
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

load("data/CARE_SS2_pData_MPstates.RData", verbose = T)
pDataSS2 <- pData

pData10x <- readRDS("data/malignant_meta_data_2025_01_08.RDS")

mDataSS2 <- table(pDataSS2$ID, pDataSS2$State)
mDataSS2 <- as.matrix(mDataSS2)
mDataSS2 <- t(apply(mDataSS2, 1, function(x) (x/sum(x))*100))

mData10x <- table(pData10x$ID, pData10x$State)
mData10x <- as.matrix(mData10x)
mData10x <- t(apply(mData10x, 1, function(x) (x/sum(x))*100))

cs <- intersect(rownames(mDataSS2), rownames(mData10x))
mDataSS2 <- mDataSS2[cs, ]
mData10x <- mData10x[cs, ]

ccor <- cor(mData10x, mDataSS2)
ccor <- ccor[rownames(ccor) != "Unresolved", colnames(ccor) != "Unresolved"]

par(mar=c(8, 5, 8, 2.1))
corrplot(ccor, method="color", is.corr=T,
         pch.cex = 3, pch.col = "white",
         col=colorRampPalette(c("blue4", "white", "firebrick"))(75),
         tl.srt = 45, tl.col = "black", tl.cex = 1, cl.pos = "r", mar = c(0, 0, 0, 5))
