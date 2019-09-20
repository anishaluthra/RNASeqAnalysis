#' Function "plotHeatmap"
#'
#' This function plots heat maps using given input matrix and group conditions
#' @importFrom Biobase exprs
#' @importFrom genefilter varFilter
#' @importFrom CancerSubtypes FSbyMAD
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @import ComplexHeatmap
#' @param geneExp a matrix, rows correspond to genes and columns correspond to samples if it is a gene expression matrix.
#' @param group a data frame with each column a different column annotation which needs to be added to the heatmap, rownames correspond to sample ID.
#' @param myPalette a named list with color palettes for each annotation, each element in the list should be a named vector.
#' @param legendName name of the legend, default is "Gene expression".
#' @param method a character type variable indicating the clustering method. Default is "ward.D2".
#' @param distance a character type variable indicating the distance used for clustering. Default is "spearman".
#' @param rowName a logical variable indicating whether to show the row names. Default is FALSE.
#' @param s1 a numerical variable indicating the fontsize of legend title. Default is 15.
#' @param s2 a numerical variable indicating the fontsize of row names. Default is 8.
#' @param s3 a numerical variable indicating the fontsize of column names. Default is 12.
#' @param s4 a numerical variable indicating the height and width of grid in defined unit of column legend. Default is 6.
#' @param s5 a numerical variable indicating the width of the dendrogram for columns and rows. Default is 12.
#' @param s6 a numerical variable indicating the width of space for rowname and colname annotations. Default is 100.
#' @param H a numerical variable indicating the height of saved figure. Default is 15.
#' @param W a numerical variable indicating the width of save figure. Default is 15.
#' @param Date a Date object obtained from Sys.Date
#' @return a "HeatmapList" class object.
#' @export

plotHeatmap <- function(geneExp, group, myPalette, legendName="Gene expression", method="ward.D2", distance="spearman",
                        rowName=FALSE, s1=15, s2=8, s3=12, s4=6, s5=12, s6=100, H, W, title="Heatmap", Date) {
  # construct color palette
  #myPalette <- list(c("#EE82EE", "#00FFFF"), c("#7D26CD", "#00CD00"), c("#FEB24C", "#FC4E2A", "#BD0026"),
  #                  c("#00BFFF", "#1E90FF", "#00008B"), c("#90EE90", "#00CD00", "#006400", "#EEAEEE"),
  #                  c("#FF69B4", "#FF1493", "#C51B7D", "#DAA520", "#999999"), c("#9932CC", "#FF8C00"))
  #names(myPalette) <- colnames(group)
  #for (i in 1:ncol(group)) {
  #  names(myPalette[[i]]) <- levels(group[, i])
  #}

  # construct column annotations
  legendCol <- vector(mode = "list")
  for (i in 1:ncol(group)) {
    legendCol[[i]] <- list(title = colnames(group)[i], at = levels(group[, colnames(group)[i]]), labels = levels(group[, colnames(group)[i]]),
                           title_gp = gpar(fontsize = s1, fontface = "bold"), labels_gp = gpar(fontsize = s1),
                           grid_height = unit(s4, "mm"), grid_width = unit(s4, "mm"))
  }
  names(legendCol) <- colnames(group)
  colAnno <- HeatmapAnnotation(df = group,
                               show_annotation_name = NULL, annotation_legend_param = legendCol,
                               col = myPalette, show_legend = FALSE)

  # plot heatmap for whole group
  ht <- Heatmap(geneExp, cluster_row_slices = FALSE, cluster_column_slices = FALSE, row_labels = row.names(geneExp),
                name = legendName, cluster_rows = TRUE, cluster_columns = TRUE, show_row_dend = TRUE, show_column_dend = TRUE,
                show_row_names = rowName, clustering_distance_rows = distance, clustering_distance_columns = distance,
                clustering_method_rows = method, clustering_method_columns = method, row_names_gp = gpar(fontsize = s2),
                column_names_rot = -45, column_names_gp = gpar(fontsize = s3), top_annotation = colAnno, row_title = NULL,
                column_title = NULL, heatmap_legend_param = list(ncol = 1, title_gp = gpar(fontsize = s1, fontface = "bold"),
                                                                 labels_gp = gpar(fonrsize = s1), grid_height = unit(s4, "mm"),
                                                                 grid_width = unit(s4, "mm")),
                show_heatmap_legend = FALSE, row_dend_width = unit(s5, "mm"), column_dend_height = unit(s5, "mm"),
                column_names_max_height = unit(s6, "mm"), row_names_max_width = unit(s6, "mm"), use_raster = FALSE)

  colLegend <- lapply(colAnno@anno_list, function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE,
                                                                             title_gp = gpar(fontsize = s1, fontface = "bold"),
                                                                             labels_gp = gpar(fontsize = s1), grid_height = unit(s4, "mm"),
                                                                             grid_width = unit(s4, "mm")))
  hmLegend <- color_mapping_legend(ht@matrix_color_mapping, plot = FALSE, title_gp = gpar(fontsize = s1, fontface = "bold"),
                                   labels_gp = gpar(fontsize = s1), grid_height = unit(s4, "mm"), grid_width = unit(s4, "mm"))
  ht <- draw(ht, heatmap_legend_list = c(colLegend, hmLegend), merge_legends = TRUE)
  g <- arrangeGrob(grid.grabExpr(draw(ht)), nrow = 1,ncol = 1)
  ggsave(g, filename=paste0("Heatmap_", title, "_", Date, ".pdf"), height = H, width = W, dpi = "retina", limitsize = FALSE)

  # return results
  return(ht)
}
