#' Function "NESplot"
#'
#' This function generates an normalized enrichment score plot.
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @param gseaRes the result from GSEA
#' @param title a string character indicating the name of geneset
#' @param type a character indicating the group condition
#' @param H height for ES plot
#' @param W width for ES plot
#' @param Date an Date object
#' @return a ggplot object of the nes plot
#' @export

NESplot <- function(gseaRes, pval, title, type, H, W, Date) {
  # plot the normalized enrichment scores
  # color the bar indicating whether or not the pathway was significant
  gseaResI <- gseaRes
  gseaResI$Group <- rep(paste0("padj<", pval), nrow(gseaResI))
  gseaResI$Group <- sapply(1:nrow(gseaResI), FUN=function(x) {
    ifelse(gseaResI$padj[x] <= pval, paste0("padj<", pval), paste0("padj>=", pval))
  }) %>% as.factor()
  myPalette <- c("#00BFC4", "#F8766D")
  names(myPalette) <- c(paste0("padj<", pval), paste0("padj>=", pval))
  myBreaks <- c(paste0("padj<", pval), paste0("padj>=", pval))
  
  nes <- ggplot(gseaResI, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=Group)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(title, " geneset NES from GSEA (", type, ")")) +
    theme_minimal() +
    scale_fill_manual(breaks=myBreaks, values=myPalette)
  
  g <- arrangeGrob(nes, nrow = 1)
  ggsave(g, filename=paste0("NESplot_", title, "_", type, "_", Date, ".pdf"), height=H, width=W, limitsize=FALSE)
  
  return(nes)
}