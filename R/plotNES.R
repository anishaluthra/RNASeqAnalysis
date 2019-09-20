#' Function "plotNES"
#'
#' This function plots a comparison pathway barplot for selected cohorts and chosen gene sets.
#' @importFrom magrittr %>%
#' @importFrom gridExtra arrangeGrob
#' @importFrom dplyr left_join
#' @import ggplot2
#' @param gseaRes a list of gsea results for different cohorts
#' @param N a numeric variable indicating the number of cohorts being used
#' @param gs a character variable indicating the name of genesets used
#' @param fdr a numeric variable indicating the FDR cutoff value
#' @param s a numeric variable representing the size of labels used in the plot
#' @param v a numeric variable representing the value of vjust, vertical position of the labels
#' @param h a numeric variable representing the value of hjust, horizontal position of the labels
#' @param W a numeric variable representing the width of the figure
#' @param H a numeric variable representing the height of the figure
#' @param Date a character variable returned from Sys.Date()
#' @return a list of ggplot object of the comparison plot of direction-consistent pathways, up and down-regulated pathways
#' @export

plotNES <- function(gseaRes, N, gs, fdr, s, v, h, W, H, Date) {
  # construct a data frame for plotting
  df <- Reduce(rbind, gseaRes)[, c("pathway", "padj")]
  df$`-log10(padj)` <- -log10(df$padj)
  df$Cohort <- rep(names(gseaRes), each = nrow(gseaRes[[1]]))

  # figure out the up-regulated and down-regulated pathways in all cohorts
  list1 <- lapply(X=1:length(gseaRes), FUN=function(x) {
    gseaRes[[x]]$pathway[gseaRes[[x]]$NES > 0]
  })
  list2 <- lapply(X=1:length(gseaRes), FUN=function(x) {
    gseaRes[[x]]$pathway[gseaRes[[x]]$NES < 0]
  })
  upPath <- Reduce(intersect, list1)
  dnPath <- Reduce(intersect, list2)

  # filter out those pathways not in the upPath or dnPath
  df1 <- data.frame(pathway = rep(c(upPath, dnPath), N), Cohort = rep(names(gseaRes), each = length(c(upPath, dnPath))))
  df1$pathway <- as.character(df1$pathway)
  df1$Cohort <- as.character(df1$Cohort)
  df <- left_join(df1, df)
  df$Group <- rep(c(rep("UP", length(upPath)), rep("DOWN", length(dnPath))), length(gseaRes))
  df$`-log10(padj)`[which(df$Group == "DOWN")] <- df$`-log10(padj)`[which(df$Group == "DOWN")] * (-1)
  df$Condition <- rep(paste0("padj < ", fdr), nrow(df))
  df$Condition <- sapply(1:nrow(df), FUN=function(x) {
    ifelse(df$padj[x] < fdr, paste0("padj < ", fdr), paste0("padj >= ", fdr))
  }) %>% as.factor()
  df$Label <- rep(" ", nrow(df))
  df$Label <- sapply(1:nrow(df), FUN=function(x) {
    ifelse(df$padj[x] < fdr, "**", " ")
  }) %>% as.factor()
  #df2 <- df[df$Cohort == names(gseaRes)[1], ]
  #df2 <- df2[order(df2$`-log10(padj)`, decreasing = TRUE), ]
  #row.names(df2) <- df2$pathway
  #df3 <- df
  #for (i in 1:N) {
  #  df4 <- df[df$Cohort == names(gseaRes)[i], ]
  #  row.names(df4) <- df4$pathway
  #  df3[(((i-1)*nrow(df2)+1):(i*nrow(df2))), ] <- df4[row.names(df2), ]
  #  row.names(df3) <- NULL
  #}

  myBreaksI <- c(paste0("padj < ", fdr), paste0("padj >= ", fdr))
  myBreaksII <- names(gseaRes)
  #myPalette <- list("#00BFC4", "#F8766D")
  myPaletteI <- c("gold", "#C0C0C0")
  myPaletteII <- as.character(brewer.pal(N, "Set1"))
  names(myPaletteI) <- myBreaksI
  names(myPaletteII) <- myBreaksII

  # plot the normalized enrichment scores
  # color the bar indicating whether or not the pathway was significant
  nes <- ggplot(df, aes(reorder(pathway, `-log10(padj)`), `-log10(padj)`, colour=Condition, width=v)) +
    geom_col(aes(fill=Cohort, colour=Condition), position=position_dodge2(preserve="single", width=1), alpha=0.7) +
    geom_text(aes(fill=Cohort, label=Label), size=s, position=position_dodge(width=0.8), vjust=v, hjust=h) +
    coord_flip() +
    labs(x="Pathway", y="-log10(padj)",
         title=paste0("GSEA Result Comparison for ", gs, " gene sets")) +
    theme_minimal() +
    scale_fill_manual(breaks=myBreaksII, values=myPaletteII) +
    scale_color_manual(breaks=myBreaksI, values=myPaletteI)
  nes <- nes + theme(legend.title = element_blank())

  g.nes <- arrangeGrob(nes, nrow = 1, ncol = 1)
  ggsave(g.nes, filename = paste0("NEScomparison_", gs, "_", Date, ".pdf"), width = W, height = H, limitsize = FALSE)

  # return results
  return(list("nes" = nes, "df" = df, "upPath" = upPath, "dnPath" = dnPath))
}
