#################################
## Title: Figure 2 panel e in Nomura et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Association between the pathway-based (rows) and gene-based (columns) MPs
#################################

# Necessary packages
library(tidyverse)
library(Matrix)
library(psych)
library(ComplexHeatmap)
library(colorspace)
library(circlize)

NES_MP <- readRDS("data/MP_scores.RDS")
tmp <- NES_MP[, 6:ncol(NES_MP)] %>%
  as.matrix
rownames(tmp) <- pull(NES_MP, CellID)
NES_MP <- tmp
rm(tmp)

NES_PMP <- readRDS("data/aggregate_pathway_scores.RDS")
cc <- intersect(rownames(NES_MP), rownames(NES_PMP))
NES_MP <- NES_MP[cc, ]
NES_PMP <- NES_PMP[cc, ]
rm(cc)

tmp <- corr.test(NES_MP, NES_PMP, method = "pearson", adjust = "none", ci = FALSE)
ccor <- tmp$r
ccorP <- tmp$p

NES_MP <- apply(NES_MP, 2, function(x) (x-mean(x))/sd(x))
ddist <- dist(t(NES_MP))
MPhc <- hclust(ddist, method = "ward.D2")
MPct <- cutree(MPhc, k = 7)

NES_PMP <- apply(NES_PMP, 2, function(x) (x-mean(x))/sd(x))
ddist <- dist(t(NES_PMP))
PMPhc <- hclust(ddist, method = "ward.D2")
PMPct <- cutree(PMPhc, k = 7)

ccolPMP <- c("forestgreen", "cyan2", "blue", "red", "orange", "magenta", "firebrick")
ccolPMP <- ccolPMP[PMPct]
names(ccolPMP) <- colnames(ccor)

ccolMP <- c("forestgreen", "blue", "cyan2", "orange", "firebrick", "red", "magenta")
ccolMP <- ccolMP[MPct]
names(ccolMP) <- rownames(ccor)

rownames(ccor) <- sapply(strsplit(rownames(ccor), "_"), function(x) paste(x[1:2], collapse = "_"))

ta <- HeatmapAnnotation(subtype = ccolMP,
                        col = list(subtype = c("firebrick" = "firebrick", "red" = "red", "forestgreen" = "forestgreen",
                                               "blue" = "blue", "cyan2" = "cyan2", "magenta" = "magenta", "orange" = "orange")
                        ),
                        show_legend = F, show_annotation_name = F
)
la <- rowAnnotation(groups = ccolPMP,
                    col = list(groups = c("firebrick" = "firebrick", "red" = "red", "forestgreen" = "forestgreen",
                                          "blue" = "blue", "cyan2" = "cyan2", "magenta" = "magenta", "orange" = "orange")
                    ),
                    show_legend = F, show_annotation_name = F
)

col_fun <- colorRamp2(seq(-.5, .5, length.out = 75), colorRampPalette(c("blue", "lightblue", "white", "orange", "red"))(75))

Heatmap(t(ccor), name = "mat",
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE, clustering_distance_rows = "euclidean", 
        clustering_method_rows = "ward.D2", cluster_columns = as.dendrogram(MPhc),
        cluster_rows = rev(as.dendrogram(PMPhc)),
        col = col_fun,
        top_annotation = ta, show_row_dend = T,
        left_annotation = la,
        show_column_names = T, show_row_names = T,
        row_dend_width = unit(20, "mm"),
        column_dend_height = unit(20, "mm"),
        width = unit(100, "mm"), height = unit(100, "mm")
)


