#' Function "plotVolcanoII"
#'
#' This function generates volcano plot
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @param fc a numeric variable representing the log2FoldChange cutoff
#' @param pval a numeric variable representing the adjusted p value cutoff
#' @param res a data frame containing the differential gene expression analysis result(log2FoldChange and padj)
#' @param attr1 a character indicating column name of log2FoldChange
#' @param attr2 a character indicating column name of adjusted p value
#' @param mytitle a string character indicating the title of the plot
#' @param group a string character indicating the group condition
#' @param UPgenes "up-regulated" pathway in the given geneset
#' @param DWgenes "down-regulated" pathway in the given geneset 
#' @param xmin the minimum of x axis
#' @param xmax the maximum of x axis
#' @param ymin the minimum of y axis
#' @param ymax the maximum of y axis
#' @param W a numeric variable indicating the width of the plot
#' @param Date a Date object obtained from Sys.Date
#' @return a list consisting of normalized gene expression matrix.
#' @export

plotVolcanoII <- function(fc, pval, res, attr1, attr2, mytitle, group, UPgenes, DWgenes, xmin, xmax, ymin, ymax, W, Date) {
  # volcano plot
  ns <- setdiff(rownames(res), c(UPgenes, DWgenes))
  res <- res[c(ns, UPgenes, DWgenes), ]
  myCol <- rep('grey70', nrow(res))
  names(myCol) <- rep('NS', nrow(res))
  myCol[which(rownames(res) %in% UPgenes)] <- "red2"
  names(myCol)[which(rownames(res) %in% UPgenes)] <- "UP"
  myCol[which(rownames(res) %in% DWgenes)] <- "royalblue"
  names(myCol)[which(rownames(res) %in% DWgenes)] <- "DOWN"

  vp <- EnhancedVolcano(res[, c(attr1, attr2)], lab = rownames(res),
                        x = attr1, y = attr2,
                        xlim = c(xmin, xmax), ylim = c(ymin, ymax),
                        title = mytitle,
                        pCutoff = pval, FCcutoff = fc, cutoffLineType = 'blank',
                        vline = c(fc, -fc), vlineType = 'solid',
                        hline = pval, hlineType = 'solid',
                        #hlineCol = c('orange', 'cyan2'),
                        transcriptPointSize = 1.2, transcriptLabSize = 4.0,
                        xlab = "log2FoldChange", ylab = "-log10(padj)",
                        subtitle = "", colAlpha = 0.5, legendPosition = "right",
                        #col=c("grey30", "royalblue", "red2", "purple"),
                        colCustom = myCol,
                        selectLab = "",
                        legendVisible = FALSE,
                        #legendLabSize = 0, legendIconSize = 2.0,
                        gridlines.major = FALSE, gridlines.minor = FALSE)

  # save volcano plots
  g.vp <- arrangeGrob(vp, nrow=1, ncol=1)
  ggsave(g.vp, filename=paste0("VolcanoPlot_", group, "_", Date, ".png"), width = W, dpi = 300)

  # return result
  return(vp)
}
