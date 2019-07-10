#' Function "ssgseaResult"
#'
#' This function stores the ssGSEA results and plots
#' @importFrom data.table fwrite as.data.table
#' @importFrom ggplot2 ggsave
#' @importFrom gridExtra arrangeGrob
#' @param ssgseaRes an ExpressionSet object of ssGSEA result
#' @param ssgsea result containing all ssGSEA plots object
#' @param gs a string character indicating the name of gene set
#' @param group a string character indicating the group condition
#' @param condition a character vector indicating the group condition
#' @param Date a Date object obtained from Sys.Date
#' @export

ssgseaResult <- function(ssgseaRes, ssgsea, gs, group, condition, Date) {
  fwrite(as.data.table(exprs(ssgseaRes), keep.rownames = TRUE),
         file = paste0("ssgseaResult_", group, "_", gs, "_", Date, ".csv"),
         sep = ",", quote = "auto", row.names = FALSE)

  gI <- arrangeGrob(ssgsea$hmI.1$gtable, ssgsea$hmI.2$gtable, ncol = 2, nrow = 1)
  ggsave(gI, filename = paste0("HeatmapSSGSEA_", group, "_", gs, "_", Date, ".pdf"), width = 60,
         height = 135, limitsize = FALSE)
  gII <- arrangeGrob(ssgsea$hmII.1$gtable, ssgsea$hmII.2$gtable, ncol = 2, nrow = 1)
  ggsave(gII, filename = paste0("HeatmapSSGSEA(", condition[1], ")_", group, "_", gs, "_", Date, ".pdf"), width = 40,
         height = 135, limitsize = FALSE)
  gIII <- arrangeGrob(ssgsea$hmIII.1$gtable, ssgsea$hmIII.2$gtable, ncol = 2, nrow = 1)
  ggsave(gIII, filename = paste0("HeatmapSSGSEA(", condition[2], ")_", group, "_", gs, "_", Date, ".pdf"), width = 50,
         height = 135, limitsize = FALSE)
  gIV <- arrangeGrob(ssgsea$bp, ncol = 1, nrow = 1)
  ggsave(gIV, filename = paste0("BoxplotSSGSEA_", group, "_", gs, "_", Date, ".pdf"), width = 5,
         dpi = 300, limitsize = FALSE)
  gV <- arrangeGrob(ssgsea$bp1, ncol = 1, nrow = 1)
  ggsave(gV, filename = paste0("BoxplotSSGSEA(pathway)_", group, "_", gs, "_", Date, ".pdf"), width = 200,
         height = 30, limitsize = FALSE)
  gVI <- arrangeGrob(ssgsea$bp2, ncol = 1, nrow = 1)
  ggsave(gVI, filename = paste0("BoxplotSSGSEA(significant pathway)_", group, "_", gs, "_", Date, ".pdf"), width = 50,
         height = 20, limitsize = FALSE)
}
