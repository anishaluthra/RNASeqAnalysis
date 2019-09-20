#' Function "plotHeatmap_POI"
#'
#' This function generates heatmap for the genes from selected pathways(pathway of interest)
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @param condition a vector of characters representing all groups of samples, default is c("mutant", "wildtype").
#' @param geneES an ExpressionSet object representing the gene expression matrix, at least two columns in pData(geneES) is required,
#'               "SampledID" indicating samples and "Group" indicating which condition the samples are in.
#' @param group a data frame with each column a different column annotation which needs to be added to the heatmap, rownames correspond to sample ID.
#' @param palette1 a named list with color palettes for each annotation, each element in the list should be a named vector.
#' @param palette2 a character type vector representing the colors used for distinguishing "POI", default is
#'                 c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A").
#' @param groupName a character variable indicating the POI name, default is "Pathway".
#' @param legendName a character variable indicating the legend name, default is "Gene expression".
#' @param POI a character type vector consisting of names of selected pathways.
#' @param GOI a list consisting of gene names in the pathways of "POI"(can be all or partial genes).
#' @param rowName a boolean variable indicating whether to show the row names. Default is TRUE.
#' @param distance a character variable indicating which distance measure to be used in clustering, either one of
#'                 "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall",
#'                 or a dist object. Default is "spearman".
#' @param method a character variable indicating the clustering method used, pass to hclust. Default is "ward.D2".
#' @param title a character variable containing information for filename.
#' @param s1 a numerical variable indicating the fontsize of legend title. Default is 15.
#' @param s2 a numerical variable indicating the fontsize of row names. Default is 4.
#' @param s3 a numerical variable indicating the fontsize of column names. Default is 11.
#' @param s4 a numerical variable indicating the height and width of grid in defined unit of column legend. Default is 6.
#' @param s5 a numerical variable indicating the width of the dendrogram for columns and rows. Default is 12.
#' @param H a numerical variable indicating the height of saved figure. Default is 11.
#' @param W a numerical variable indicating the width of save figure. Default is 28.
#' @param Date a Date object obtained from Sys.Date
#' @return a "HeatmapList" class object.
#' @export


plotHeatmap_POI <- function(condition=c("mutant", "wildtype"), geneES, group, palette1,
                            palette2=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"),
                            groupName="Pathway", legendName="Gene expression", POI, GOI, rowName=TRUE, distance="spearman", method="ward.D2", title,
                            s1=15, s2=4, s3=11, s4=6, s5=12, H=11, W=28, Date) {
  # construct new gene expression matrix with genes for each pathway in POI separately
  len <- length(POI)
  geneExp <- as.data.frame(exprs(geneES)[GOI[[1]], ])
  geneExp$SampleID <- row.names(geneExp)
  row.names(geneExp) <- NULL
  if (len >= 2) {
    for (i in 2:length(POI)) {
      exp_poi <- as.data.frame(exprs(geneES[GOI[[i]], ]))
      exp_poi$SampleID <- row.names(exp_poi)
      row.names(exp_poi) <- NULL
      geneExp <- rbind(geneExp, exp_poi)
    }
  }
  geneExp_new <- as.matrix(geneExp[, 1:ncol(exprs(geneES))])
  row.names(geneExp_new) <- as.character(1:nrow(geneExp_new))

  # construct row and column split groups(gaps among different pathways)
  len_poi <- unlist(Map(length, GOI))
  rowSplit <- factor(rep(POI, len_poi), levels = POI)
  colSplit <- factor(pData(geneES)$Group, levels = condition)

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
                               col = palette1, show_legend = FALSE, gap = unit(3, "points"))


  #myColor <- list(palette1)
  #names(myColor) <- conditionName
  #names(myColor[[conditionName]]) <- condition
  #df1 <- as.data.frame(as.character(pData(geneES)$Group))
  #colnames(df1) <- conditionName
  #legendCol <- list(list(title = conditionName, at = condition, labels = condition, title_gp = gpar(fontsize = s1, fontface = "bold")))
  #names(legendCol) <- conditionName
  #colAnno <- HeatmapAnnotation(df = df1, show_annotation_name = NULL, annotation_legend_param = legendCol,
  #                             col = myColor, show_legend = FALSE, gap = unit(3, "points"))

  # construct row annotations
  myPalette <- list(palette2)
  names(myPalette) <- groupName
  names(myPalette[[groupName]]) <- POI
  df2 <- as.data.frame(as.character(rowSplit))
  colnames(df2) <- groupName
  legendRow <- list(list(title = groupName, at = POI, labels = POI, title_gp = gpar(fontsize = s1, fontface = "bold")))
  names(legendRow) <- groupName
  rowAnno <- rowAnnotation(df = df2, annotation_legend_param = legendRow,
                           col = myPalette, show_legend = FALSE, gap = unit(3, "points"))

  # plot heatmap for whole group
  ht <- Heatmap(geneExp_new, row_split = rowSplit, column_split = colSplit, cluster_row_slices = FALSE,
                cluster_column_slices = FALSE, row_labels = geneExp$SampleID, name = legendName, cluster_rows = TRUE,
                cluster_columns = TRUE, show_row_dend = TRUE, show_column_dend = TRUE, show_row_names = rowName,
                clustering_distance_rows = distance, clustering_distance_columns = distance, clustering_method_rows = method,
                clustering_method_columns = method, row_names_gp = gpar(fontsize = s2), column_names_rot = -45,
                column_names_gp = gpar(fontsize = s3), row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
                top_annotation = colAnno, left_annotation = rowAnno, row_title = NULL, column_title = NULL,
                show_heatmap_legend = FALSE, row_dend_width = unit(s5, "mm"), column_dend_height = unit(s5, "mm"),
                heatmap_legend_param = list(ncol = 1, title_gp = gpar(fontsize = s1, fontface = "bold"),
                                            labels_gp = gpar(fonrsize = s1), grid_height = unit(s4, "mm"),
                                            grid_width = unit(s4, "mm")), use_raster = FALSE)

  colLegend <- lapply(colAnno@anno_list, function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE,
                                                                             title_gp = gpar(fontsize = s1, fontface = "bold"),
                                                                             labels_gp = gpar(fontsize = s1), grid_height = unit(s4, "mm"),
                                                                             grid_width = unit(s4, "mm")))
  rowLegend <- lapply(rowAnno@anno_list[groupName], function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE,
                                                                                        title_gp = gpar(fontsize = s1, fontface = "bold"),
                                                                                        labels_gp = gpar(fontsize = s1), grid_height = unit(s4, "mm"),
                                                                                        grid_width = unit(s4, "mm")))
  hmLegend <- color_mapping_legend(ht@matrix_color_mapping, plot = FALSE, title_gp = gpar(fontsize = s1, fontface = "bold"),
                                   labels_gp = gpar(fontsize = s1), grid_height = unit(s4, "mm"), grid_width = unit(s4, "mm"))
  ht <- draw(ht, heatmap_legend_list = c(list(hmLegend), colLegend, rowLegend), merge_legends = TRUE,
             legend_title_gp = gpar(fontsize = s1, fontface = "bold"), legend_labels_gp = gpar(fontsize = s1))

  #colLegend <- lapply(colAnno@anno_list[conditionName], function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
  g <- arrangeGrob(grid.grabExpr(draw(ht)), nrow = 1,ncol = 1)
  ggsave(g, filename=paste0("HeatmapPOI_", title, "_", Date, ".pdf"), height = H, width = W, limitsize = FALSE)

  # return results
  return(ht)
}
