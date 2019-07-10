#' Function "ssgseaPlot"
#'
#' This function performs single sample GSEA
#' @importFrom stats wilcox.test
#' @importFrom pheatmap pheatmap
#' @import ggplot2
#' @importFrom Biobase exprs pData
#' @param ssgsea an ExpressionSet object containing the result of ssGSEA
#' @param condition a string character representing the group condition
#' @param n1 a numeric variable indicating the number of samples for one of the group condition
#' @param n2 a numeric variable indicating the number of samples for the other group condition
#' @param data_info data frame containing clinical information
#' @param p a numeric variable indicating the p value cutoff of wilcoxon test
#' @param attr1 a string character indicating the column names in clinical file for SampleID
#' @param attr2 a string character indicating the column containing group information in clinical data
#' @param conditionI a character vector of alternative group condition
#' @return a list of plots and p values of wilcoxon test
#' @export

ssgseaPlot <- function(ssgsea, condition, n1, n2, data_info, p, attr1, attr2, conditionI) {
  # rank test for each pathway
  rtest <- sapply(X=1:nrow(exprs(ssgsea)), FUN=function(k) {
    round(wilcox.test(x = exprs(ssgsea)[k, 1:n1], y = exprs(ssgsea)[k, (n1+1):(n1+n2)],
                      alternative = "two.sided", mu = 0, conf.level = 0.95)$p.value, digits=4)})

  rtestI <- rtest
  rtestI <- sapply(X=1:length(rtest), FUN=function(k) {
    rtestI[k] <- ifelse(as.numeric(rtest[k]) == 0, "< 0.0001", rtest[k])
  })

  # plot heatmap using whole group
  group <- data.frame(group=pData(ssgsea)$Group, condition=rep(0, ncol(exprs(ssgsea))))
  row.names(group) <- colnames(exprs(ssgsea))
  rnInfo <- data_info[which(data_info$SampleID %in% rownames(group)), c(attr1, attr2)]
  group$condition <- lapply(1:nrow(group), FUN=function(x) {
    group$condition[x] <- rnInfo[, attr2][which(rnInfo[, attr1]==rownames(group)[x])]
  })
  for (i in 1:nrow(group)) {
    ifelse(group$condition[i] == 1,
           group$condition[i] <- conditionI[1],
           group$condition[i] <- conditionI[2])
  }
  group$condition <- as.character(group$condition)
  group <- group[, c(2, 1)]

  # construct annotation color palette
  myColor <- c("violet", "cyan", "purple3", "green3")
  ann_color <- list(group = c(myColor[1], myColor[2]),
                    condition = c(myColor[3], myColor[4]))
  names(ann_color$group) <- condition
  names(ann_color$condition) <- conditionI

  hmI.1 <- pheatmap(mat=ssgsea, angle_col=45,
                    main="", annotation_col=group, annotation_colors=ann_color,
                    cluster_rows=TRUE, cellwidth=14, cellheight=16)
  hmI.2 <- pheatmap(mat=ssgsea, scale = "row", angle_col=45,
                    main="", annotation_col=group, annotation_colors=ann_color,
                    cluster_rows=TRUE, cellwidth=14, cellheight=16)

  # use only mutation group
  hmII.1 <- pheatmap(mat=ssgsea[, ssgsea$Group==condition[1]],
                     main="", annotation_col=group, angle_col=45, annotation_colors=ann_color,
                     cluster_rows=TRUE, cellwidth=14, cellheight=16)
  hmII.2 <- pheatmap(mat=ssgsea[, ssgsea$Group==condition[1]], scale = "row",
                     main="", annotation_col=group, angle_col=45, annotation_colors=ann_color,
                     cluster_rows=TRUE, cellwidth=14, cellheight=16)

  # use only wildtype group
  hmIII.1 <- pheatmap(mat=ssgsea[, ssgsea$Group==condition[2]],
                      main="", annotation_col=group, angle_col=45, annotation_colors=ann_color,
                      cluster_rows=TRUE, cellwidth=14, cellheight=16)
  hmIII.2 <- pheatmap(mat=ssgsea[, ssgsea$Group==condition[2]], scale = "row",
                      main="", annotation_col=group, angle_col=45, annotation_colors=ann_color,
                      cluster_rows=TRUE, cellwidth=14, cellheight=16)

  # boxplot mutation versus wildtype for whole geneset
  df <- data.frame(ssgseaScore = as.vector(exprs(ssgsea)),
                   Group = c(rep(condition[1], nrow(exprs(ssgsea))*n1),
                             rep(condition[2], nrow(exprs(ssgsea))*n2)))
  bp <- ggplot(df, aes(x=df$Group, y=df$ssgseaScore, fill=df$Group)) +
    geom_boxplot(alpha=0.6) +
    theme_gray() +
    ggtitle("") +
    labs(x="group", y="ssgsea enrichment score") +
    #scale_y_continuous(breaks = seq(0, 4000, 250), limits = c(0, 4000)) +
    scale_x_discrete(breaks = c(condition[1], condition[2]),
                     labels = c(condition[1], condition[2]))
  #facet_wrap(~Ensembl_id, nrow=2, scales="free") +
  bp <- bp + theme(legend.text=element_text(size=14, vjust=1, margin=margin(6, 6, 6, 6)),
                   axis.text=element_text(size=14))

  # boxplot per pathway
  myPalette <- rep("black", nrow(exprs(ssgsea)))
  for (i in 1:length(myPalette)) {
    if (rtest[i] <= p) {
      myPalette[i] <- "red"
    }
  }
  df1 <- data.frame(ssgseaScore = as.vector(t(exprs(ssgsea))),
                    Pathway = rep(rownames(exprs(ssgsea)), each = ncol(exprs(ssgsea))),
                    Group = rep(c(rep(condition[1], n1), rep(condition[2], n2)),
                                nrow(exprs(ssgsea))))
  bp1 <- ggplot(df1, aes(x=df1$Pathway, y=df1$ssgseaScore, fill=df1$Group)) +
    geom_boxplot(alpha=0.6) +
    theme_gray() +
    ggtitle("") +
    labs(x="pathway", y="ssgsea enrichment score") +
    #scale_y_continuous(breaks = seq(-4000, 8000, 1000), limits = c(-4000, 8000)) +
    scale_x_discrete(breaks = rownames(exprs(ssgsea)),
                     labels = paste0(rownames(exprs(ssgsea)), "(p: ", rtestI, ")"))
  #facet_wrap(~Ensembl_id, nrow=2, scales="free") +
  bp1 <- bp1 + theme(axis.text.x = element_text(angle = 270, hjust = 0, color = myPalette),
                     #legend.title=element_blank(),
                     legend.text=element_text(size=14, vjust=1, margin=margin(6, 6, 6, 6)),
                     axis.text=element_text(size=14))

  # only plot significant pathway boxplots(either p value < 0.05)
  indices <- which(rtest <= 0.05)
  myPaletteI <- rep("black", length(indices))
  ssgseaSig <- exprs(ssgsea)[indices, ]
  ssgseaSigI <- data.frame(Pathway = rownames(ssgseaSig), Pvalue = rtestI[indices])
  ssgseaSigI <- ssgseaSigI[order(ssgseaSigI$Pvalue), ]
  ssgseaSig <- ssgseaSig[as.numeric(rownames(ssgseaSigI)), ]
  df2 <- data.frame(ssgseaScore = as.vector(t(ssgseaSig)),
                    Pathway = rep(rownames(ssgseaSig), each = ncol(ssgseaSig)),
                    Group = rep(c(rep(condition[1], n1), rep(condition[2], n2)), nrow(ssgseaSig)))
  df2$Pathway <- factor(df2$Pathway, levels = rownames(ssgseaSig),ordered = TRUE)
  bp2 <- ggplot(df2, aes(x=df2$Pathway, y=df2$ssgseaScore, fill=df2$Group)) +
    geom_boxplot(alpha=0.6) +
    theme_gray() +
    ggtitle("") +
    labs(x="pathway", y="ssgsea enrichment score") +
    #scale_y_continuous(breaks = seq(-4000, 8000, 1000), limits = c(-4000, 8000)) +
    scale_x_discrete(breaks = rownames(ssgseaSig),
                     labels = paste0(rownames(ssgseaSig), " (p: ", rtestI[indices][as.numeric(rownames(ssgseaSigI))], ")"))
  #facet_wrap(~Ensembl_id, nrow=2, scales="free") +
  bp2 <- bp2 + theme(axis.text.x = element_text(angle = 270, hjust = 0, color = myPaletteI),
                     #legend.title=element_blank(),
                     legend.text=element_text(size=14, vjust=1, margin=margin(6, 6, 6, 6)),
                     axis.text=element_text(size=14))

  # scatter plot X: ssgsea score & Y: gene expression
  #geneExp <- matrix(0, nrow = nrow(exprs(ssgsea)), ncol = ncol(exprs(ssgsea)))
  #for (i in 1:nrow(geneExp)) {
  #  for (j in 1:ncol(geneExp)) {
  #    count <- 0
  #    for (gene in geneset[[rownames(exprs(ssgsea))[i]]]) {
  #      if (gene %in% rownames(exprs(gene.es))) {
  #        geneExp[i, j] <- geneExp[i, j] + exprs(gene.es)[gene, j]
  #        count <- count + 1
  #      }
  #    }
  #    geneExp[i, j] <- geneExp[i, j] / count
  #  }
  #}
  #df2 <- df1
  #df2$geneExpression <- as.vector(t(geneExp))
  #sp <- ggplot(df2, aes(x=ssgseaScore, y=geneExpression)) +
  #  geom_point(aes(color=Pathway, shape=Group), alpha=0.6, size=2.5) +
  #  theme_classic() +
    #scale_color_manual(name="Legend", breaks=Pathway, values=col) +
    #scale_shape_manual(name="Legend", breaks=Group, values=sha) +
  #  labs(x="ssgsea score", y="gene expression")
  #coord_equal() +
  #scale_x_continuous(breaks=seq(0, 13, 1), limits=c(0, 13)) +
  #scale_y_continuous(breaks=seq(-6, 8, 2), limits=c(-6, 8))
  #sp <- sp + guides(col = guide_legend(ncol = 1))

  # return results
  return(list("hmI.1" = hmI.1, "hmI.2" = hmI.2, "hmII.1" = hmII.1,
              "hmII.2" = hmII.2, "hmIII.1" = hmIII.1, "hmIII.2" = hmIII.2, "bp" = bp,
              "bp1" = bp1, "bp2" = bp2, "rtest" = rtest, "rtestI" = rtestI))
}
