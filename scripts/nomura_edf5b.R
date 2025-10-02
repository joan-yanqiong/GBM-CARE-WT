#################################
## Title: Extended Data Figure 5 panel b in Nomura et al
## Date: 2025.02.15
## Author: Luciano Garofano
## Description: Clustering of the Pathway-based Meta-Programs scores based on their pairwise Pearson correlations
#################################

# Necessary packages
library(tidyverse)
library(Matrix)
library(psych)
library(ComplexHeatmap)
library(colorspace)
library(circlize)

NES <- readRDS("data/aggregate_pathway_scores.RDS")
tmp <- corr.test(NES, method = "pearson", adjust = "none", ci = FALSE)
ccor <- tmp$r
ccorP <- tmp$p

NES <- apply(NES, 2, function(x) (x-mean(x))/sd(x))
ddist <- dist(t(NES))
aHc <- hclust(ddist, method = "ward.D2")
aCT <- cutree(aHc, k = 7)

ccol <- c("forestgreen", "cyan2", "blue", "red", "orange", "magenta", "firebrick")
ccol <- ccol[aCT]
names(ccol) <- colnames(ccor)


ta <- HeatmapAnnotation(subtype = ccol,
                        col = list(subtype = c("firebrick" = "firebrick", "red" = "red", "forestgreen" = "forestgreen",
                                               "blue" = "blue", "cyan2" = "cyan2", "magenta" = "magenta", "orange" = "orange")
                        ),
                        show_legend = F, show_annotation_name = F
)
la <- rowAnnotation(groups = ccol,
                    col = list(groups = c("firebrick" = "firebrick", "red" = "red", "forestgreen" = "forestgreen",
                                          "blue" = "blue", "cyan2" = "cyan2", "magenta" = "magenta", "orange" = "orange")
                    ),
                    show_legend = F, show_annotation_name = F
)

col_fun <- colorRamp2(seq(-1, 1, length.out = 75), colorRampPalette(c("blue", "lightblue", "white", "orange", "red"))(75))

Heatmap(ccor, name = "mat",
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE, clustering_distance_rows = "euclidean", 
        clustering_method_rows = "ward.D2", cluster_columns = as.dendrogram(aHc),
        cluster_rows = as.dendrogram(aHc),
        col = col_fun,
        top_annotation = ta, show_row_dend = T,
        left_annotation = la,
        show_column_names = T, show_row_names = T,
        row_dend_width = unit(20, "mm"),
        column_dend_height = unit(20, "mm"),
        width = unit(100, "mm"), height = unit(100, "mm")
)