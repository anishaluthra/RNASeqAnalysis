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
#' @param palette1 a character type vector representing the colors used for distinguishing "condition", default is c("#EE82EE", "#00FFFF").
#' @param palette2 a character type vector representing the colors used for distinguishing "POI", default is
#'                 c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A").
#' @param conditionName a character variable indicating the condition name.
#' @param groupName a character variable indicating the POI name, default is "Pathway".
#' @param legendName a character variable indicating the legend name, default is "Gene expression".
#' @param POI a character type vector consisting of names of selected pathways.
#' @param GOI a list consisting of gene names in the pathways of "POI"(can be all or partial genes).
#' @param distance a character variable indicating which distance measure to be used in clustering, either one of
#'                 "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall",
#'                 or a dist object. Default is "spearman".
#' @param method a character variable indicating the clustering method used, pass to hclust. Default is "ward.D2".
#' @param title a character variable containing information for filename.
#' @param s1 a numerical variable indicating the fontsize of legend title. Default is 15.
#' @param s2 a numerical variable indicating the fontsize of row names. Default is 4.
#' @param s3 a numerical variable indicating the fontsize of column names. Default is 11.
#' @param H a numerical variable indicating the height of saved figure. Default is 11.
#' @param W a numerical variable indicating the width of save figure. Default is 28.
#' @param Date a Date object obtained from Sys.Date
#' @return a "HeatmapList" class object.
#' @export


plotHeatmap_POI <- function(condition=c("mutant", "wildtype"), geneES, palette1=c("#EE82EE", "#00FFFF"),
                            palette2=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"),
                            conditionName, groupName="Pathway", legendName="Gene expression", POI, GOI, distance="spearman", method="ward.D2", title,
                            s1=15, s2=4, s3=11, H=11, W=28, Date) {
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
  myColor <- list(palette1)
  names(myColor) <- conditionName
  names(myColor[[conditionName]]) <- condition
  df1 <- as.data.frame(as.character(pData(geneES)$Group))
  colnames(df1) <- conditionName
  legendCol <- list(list(title = conditionName, at = condition, labels = condition, title_gp = gpar(fontsize = s1, fontface = "bold")))
  names(legendCol) <- conditionName
  colAnno <- HeatmapAnnotation(df = df1, show_annotation_name = NULL, annotation_legend_param = legendCol,
                               col = myColor, show_legend = FALSE, gap = unit(3, "points"))

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
  plotHM1 <- Heatmap(geneExp_new, row_split = rowSplit, column_split = colSplit, cluster_row_slices = FALSE,
                     cluster_column_slices = FALSE, row_labels = geneExp$SampleID, name = legendName, cluster_rows = TRUE,
                     cluster_columns = TRUE, show_row_dend = TRUE, show_column_dend = TRUE, show_row_names = TRUE,
                     clustering_distance_rows = distance, clustering_distance_columns = distance, clustering_method_rows = method,
                     clustering_method_columns = method, row_names_gp = gpar(fontsize = s2), column_names_rot = -45,
                     column_names_gp = gpar(fontsize = s3), row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
                     top_annotation = colAnno, left_annotation = rowAnno, row_title = NULL, column_title = NULL,
                     heatmap_legend_param = list(ncol = 1), show_heatmap_legend = FALSE)
  colLegend <- lapply(colAnno@anno_list[conditionName], function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
  rowLegend <- lapply(rowAnno@anno_list[groupName], function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
  hmLegend <- color_mapping_legend(plotHM1@matrix_color_mapping, plot = FALSE)
  ht <- draw(plotHM1, heatmap_legend_list = c(colLegend, rowLegend, list(hmLegend)),
             legend_title_gp = gpar(fontsize = s1, fontface = "bold"), legend_labels_gp = gpar(fontsize = s1))
  g1 <- arrangeGrob(grid.grabExpr(draw(ht)), nrow = 1,ncol = 1)
  ggsave(g1, filename=paste0("HeatmapPOI_", title, "_", Date, ".pdf"), height = H, width = W, limitsize = FALSE)

  # return results
  return(ht)
}
