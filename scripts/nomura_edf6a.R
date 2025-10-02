#################################
## Title: Extended Data Figure 6 panel a in Nomura et al - Jaccard similiarity matrix of PCA signatures
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a Jaccard similiarity heatmap of PCA signatures
#################################

Heatmap(jac_m, use_raster = F, raster_by_magick = T,
        col = circlize::colorRamp2(c(0, .03, .06, .15, .2, .25), c("white", "white", "gold", "orange", "red", "darkred")),
        show_row_names = T, show_column_names = T,
        cluster_columns = hc, cluster_rows = hc,
        clustering_distance_rows = "euclidean", clustering_method_rows = "average",
        clustering_distance_columns = "euclidean", clustering_method_columns = "average",
        column_title = "Programs", column_title_gp = gpar(fontsize = 20),
        row_title = "Programs", row_title_gp = gpar(fontsize = 20), 
        row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(fontsize = 14),
        top_annotation = HeatmapAnnotation(Cluster = clusters,
                                           col = list(Cluster = setNames(brewer_pal(palette = "Spectral")(length(levels(clusters))), levels(clusters))),
                                           show_annotation_name = FALSE, #show_legend = FALSE,
                                           annotation_legend_param = list(title_gp = gpar(fontsize = 20),
                                                                          labels_gp = gpar(fontsize = 14), border = "black",
                                                                          grid_height = unit(7, "mm"))),
        right_annotation = rowAnnotation(Cluster = clusters,
                                         col = list(Cluster = setNames(brewer_pal(palette = "Spectral")(length(levels(clusters))), levels(clusters))),
                                         show_annotation_name = FALSE, show_legend = FALSE,
                                         annotation_legend_param = list(title_gp = gpar(fontsize = 20),
                                                                        labels_gp = gpar(fontsize = 14), border = "black",
                                                                        grid_height = unit(7, "mm"))),
        heatmap_legend_param = list(at = c(0, .05, .1, .15, .2, .25), title = "Jaccard",
                                    labels = c("0%", "5%", "10%", "15%", "20%", "25%"),
                                    title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 14), title_position = "topleft", border = "black",
                                    grid_height = unit(7, "mm")),
        column_title_side = "bottom",
        row_gap = unit(0,"mm"), cluster_column_slices = T, cluster_row_slices = F)
