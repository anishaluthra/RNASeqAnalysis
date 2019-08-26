#' Function "plotDensity"
#'
#' This function generates scatter plot with marginal density distributions
#' @importFrom ggpubr ggscatterhist
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @param info a data frame with required clinical information
#' @param ssgseaRes the result matrix of ssGSEA
#' @param axis1 a character variable indicating the label for x axis
#' @param axis2 a character variable indicating the label for y axis
#' @param type type of data, group information
#' @param Date a Date object obtained from Sys.Date
#' @return the density plot object
#' @export

plotDensity <- function(info, ssgseaRes, axis1 = "UPssgsea", axis2 = "DOWNssgsea", type, Date) {
  #groupI <- info$SampleID[which(info$Group == group1)]
  #groupII <- info$SampleID[which(info$Group == group2)]
  #groupIII <- info$SampleID[which(info$Group == group3)]
  df <- data.frame(SampleID = info$SampleID, Group = info$Group,
                   UPssgsea = ssgseaRes[rownames(ssgseaRes)[1], ],
                   DNssgsea = ssgseaRes[rownames(ssgseaRes)[2], ])
  #df$Group <- as.character(df$Group)
  #for (i in 1:ncol(ssgseaRes)) {
  #  if (as.character(df$SampleID[i]) %in% groupI) {
  #    df$Group[i] <- group1
  #  } else if (as.character(df$SampleID[i]) %in% groupII) {
  #    df$Group[i] <- group2
  #  } else if (as.character(df$SampleID[i]) %in% groupIII) {
  #    df$Group[i] <- group3
  #  }
  #}

  # construct color palette
  myColorI <- c("purple", "orange")
  myColorII <- c("red", "purple", "orange")
  if (length(unique(info$Group))==2) {
    Palette = myColorI
  } else {
    Palette = myColorII
  }
  names(Palette) <- as.character(unique(info$Group))

  plotDen <- ggscatterhist(df, x = "UPssgsea", y = "DNssgsea",
                           color = "Group", size = 3, alpha = 0.6,
                           palette = Palette,
                           #xlim = c(Min, Max), ylim = c(Min, Max),
                           margin.params = list(fill = "Group", color = "black", size = 0.2),
                           legend = "right", xlab = axis1, ylab = axis2)
  #geom_text(x = 0.7, y = 0.92, label = paste0("p_UP=", rtest[1])) +
  #geom_text(x = 0.85, y = 0.12, label = paste0("p_DOWN=", rtest[2]), angle = 270)

  g <- arrangeGrob(plotDen, nrow = 1, ncol = 1)
  ggsave(g, filename = paste0("DensityPlot_", type, "_", Date, ".pdf"), width = 10, dpi = 300, limitsize = FALSE)

  # return result
  return(plotDen)
}
