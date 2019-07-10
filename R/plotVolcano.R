#' Function "plotVolcano"
#'
#' This function generates volcano plot
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @param fc a numeric variable representing the log2FoldChange cutoff
#' @param pval a numeric variable representing the adjusted p value cutoff
#' @param res a data frame containing the differential gene expression analysis result(log2FoldChange and padj)
#' @param mytitle a string character indicating the title of the plot
#' @param group a string character indicating the group condition
#' @param xmin the minimum of x axis
#' @param xmax the maximum of x axis
#' @param ymin the minimum of y axis
#' @param ymax the maximum of y axis
#' @param W a numeric variable indicating the width of the plot
#' @param Date a Date object obtained from Sys.Date
#' @return a list consisting of normalized gene expression matrix.
#' @export

plotVolcano <- function(fc, pval, res, mytitle, group, xmin, xmax, ymin, ymax, W, Date) {
  # volcano plot
  myCol <- rep('grey70', nrow(res))
  names(myCol) <- rep('NS', nrow(res))
  myCol[which((res$log2FoldChange > fc) & (res$padj < pval))] <- "red2"
  names(myCol)[which((res$log2FoldChange > fc) & (res$padj < pval))] <- "UP"
  myCol[which((res$log2FoldChange < -fc) & (res$padj < pval))] <- "royalblue"
  names(myCol)[which((res$log2FoldChange < -fc) & (res$padj < pval))] <- "DOWN"

  vp <- EnhancedVolcano(res[, c("log2FoldChange", "padj")], lab = rownames(res),
                        x = "log2FoldChange", y = "padj",
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
                        #selectLab = "ENSG00000167244",
                        selectLab = "", legendVisible = FALSE,
                        #legendLabSize = 0, legendIconSize = 2.0,
                        gridlines.major = FALSE, gridlines.minor = FALSE)

  # save volcano plots
  g.vp <- arrangeGrob(vp, nrow=1, ncol=1)
  ggsave(g.vp, filename=paste0("VolcanoPlot_", group, "_", Date, ".png"), dpi = 300, width = W)

  # return result
  return(vp)
}
