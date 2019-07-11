#' Function "plotGSEA"
#'
#' This function generates variety of plots from GSEA result
#' @importFrom magrittr %>%
#' @importFrom pheatmap pheatmap
#' @import ggplot2
#' @importFrom Biobase exprs pData
#' @param data_info data frame containing clinical information
#' @param gseaRes the result from GSEA
#' @param gs geneset "CORALINE Signature" generated from our cohort
#' @param gs_name the geneset used for GSEA based on gene names
#' @param geneRank the gene ranking used for GSEA
#' @param gene.esI expression set of normalized gene counts matrix based on gene names
#' @param title a string character indicating the name of geneset
#' @param K number of top enriched/significant pathways to be plotted
#' @param type a character indicating the group condition
#' @param attr1 a character indicating the column name of sample IDs in the clinical data
#' @param attr2 a character indicating the column name of desired group in the clinical data
#' @param conditionI a character vector indicating the primary group condition
#' @param conditionII a character vector indicating the secondary group condition
#' @param N a variable indicating one of the value from attr2
#' @param H1 height for plotting significant pathways
#' @param W1 width for plotting significant pathway
#' @param H2 height for plotting enriched pathways
#' @param W2 width for plotting enriched pathways
#' @return a list of plot object
#' @export

plotGSEA <- function(data_info, gseaRes, gs, gs_name, geneRank, gene.esI, title, K, type, attr1, attr2,
                     conditionI, conditionII, N, H1, W1, H2, W2) {
  # plot the normalized enrichment scores
  # color the bar indicating whether or not the pathway was significant
  nes <- ggplot(gseaRes, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.1)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(title, " geneset NES from GSEA (", type, ")")) +
    theme_minimal()

  # plot heatmap for genes from the intersection of top K significant(least pval) pathways and "CORALINE Signature"
  # extract all genes from "CORALINE Signature"
  genes_coraline <- unlist(gs) %>% as.vector() %>% unique() %>% as.character()

  # construct group information
  group <- data.frame(group=pData(gene.esI)$Group, condition=rep(0, ncol(exprs(gene.esI))))
  row.names(group) <- colnames(exprs(gene.esI))
  rnInfo <- data_info[which(data_info[, attr1] %in% row.names(group)), c(attr1, attr2)]
  group$condition <- lapply(1:nrow(group), FUN=function(x) {
    group$condition[x] <- rnInfo[, attr2][which(rnInfo[, attr1]==row.names(group)[x])]
  })
  for (i in 1:nrow(group)) {
    ifelse(group$condition[i] == N,
           group$condition[i] <- conditionII[1],
           group$condition[i] <- conditionII[2])
  }
  group$condition <- as.character(group$condition)
  group <- group[, c(2, 1)]

  # construct annotation color palette
  myColor <- c("violet", "cyan", "purple3", "green3")
  ann_color <- list(group = c(myColor[1], myColor[2]),
                    condition = c(myColor[3], myColor[4]))
  names(ann_color$group) <- conditionI
  names(ann_color$condition) <- conditionII

  #hmSigI <- vector(mode = "list")
  #hmSigII <- vector(mode = "list")
  #for (i in 1:K) {
  #  pathway <- gseaRes[order(gseaRes$pval), ]$pathway[i]
  #  genes <- gs_name[[pathway]]
  #  if (length(genes) > 1) {
  #    hmSigI[[i]] <- pheatmap(mat=gene.esI[genes, ],
  #                            main=paste0("Genes from ", i, " Significant Pathway: ", pathway),
  #                            angle_col=45, annotation_col=group, cluster_rows=TRUE,
  #                            annotation_colors=ann_color, cellwidth=14, cellheight=14)
  #    hmSigII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
  #                             main=paste0("Genes from ", i, " Significant Pathway: ", pathway),
  #                             angle_col=45, annotation_col=group, cluster_rows=TRUE,
  #                             annotation_colors=ann_color, cellwidth=14, cellheight=14, )
  #  } else {
  #    hmSigI[[i]] <- pheatmap(mat=gene.esI[genes, ],
  #                            main=paste0("Genes from ", i, " Significant Pathway: ", pathway),
  #                            angle_col=45, annotation_col=group, cluster_rows=FALSE,
  #                            annotation_colors=ann_color, cellwidth=14, cellheight=14)
  #    hmSigII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
  #                             main=paste0("Genes from ", i, " Significant Pathway: ", pathway),
  #                             angle_col=45, annotation_col=group, cluster_rows=FALSE,
  #                             annotation_colors=ann_color, cellwidth=14, cellheight=14, )
  #  }
    ## do not cluster rows and preserve the gene order of the first heatmap
    #rowOrder <- hmSigI[[i]]$tree_row$order
    #hmSigII[[i]] <- pheatmap(mat=gene.esI[genes, ][rowOrder, ], scale="row",
    #                         main=paste0("Genes from ", i, " Significant Pathway: ", pathway),
    #                         angle_col=45, annotation_col=group, cluster_rows=FALSE,
    #                         cellwidth=14, cellheight=16, )
  #}

  # plot heatmap for the "core" genes from the intersection of top K significant pathways and "CORALINE Signature"
  hmCoreSigI <- vector(mode = "list")
  hmCoreSigII <- vector(mode = "list")
  for (i in 1:K) {
    pathway <- gseaRes[order(gseaRes$pval), ]$pathway[i]
    genes <- gseaRes[order(gseaRes$pval), ]$leadingEdgeGeneNames[[i]]
    genes <- intersect(genes, genes_coraline)
    if (length(genes) > 1) {
      hmCoreSigI[[i]] <- pheatmap(mat=gene.esI[genes, ],
                                  main=paste0("Core Genes from ", i, " Significant Pathway: ", pathway),
                                  angle_col=45, annotation_col=group, cluster_rows=TRUE,
                                  annotation_colors=ann_color, cellwidth=14, cellheight=14)
      hmCoreSigII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
                                   main=paste0("Core Genes from ", i, " Significant Pathway: ", pathway),
                                   angle_col=45, annotation_col=group, cluster_rows=TRUE,
                                   annotation_colors=ann_color, cellwidth=14, cellheight=14, )
    } else {
      hmCoreSigI[[i]] <- pheatmap(mat=gene.esI[genes, ],
                                  main=paste0("Core Genes from ", i, " Significant Pathway: ", pathway),
                                  angle_col=45, annotation_col=group, cluster_rows=FALSE,
                                  annotation_colors=ann_color, cellwidth=14, cellheight=14)
      hmCoreSigII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
                                   main=paste0("Core Genes from ", i, " Significant Pathway: ", pathway),
                                   angle_col=45, annotation_col=group, cluster_rows=FALSE,
                                   annotation_colors=ann_color, cellwidth=14, cellheight=14, )
    }
    #rowOrder <- hmCoreSigI[[i]]$tree_row$order
    #hmCoreSigII[[i]] <- pheatmap(mat=gene.esI[genes, ][rowOrder, ], scale="row",
    #                             main=paste0("Core Genes from ", i, " Significant Pathway: ", pathway),
    #                             angle_col=45, annotation_col=group, cluster_rows=FALSE,
    #                             cellwidth=14, cellheight=16, )
  }

  # plot heatmap for all the genes from the top K enriched(largest absolute NES) pathways, respecitvely
  #hmEnrI <- vector(mode = "list")
  #hmEnrII <- vector(mode = "list")
  #for (i in 1:K) {
  #  pathway <- gseaRes[order(abs(gseaRes$NES), decreasing=TRUE), ]$pathway[i]
  #  genes <- gs_name[[pathway]]
  #  if (length(genes) > 1) {
  #    hmEnrI[[i]] <- pheatmap(mat=gene.esI[genes, ],
  #                            main=paste0("Genes from ", i, " Enriched Pathway: ", pathway),
  #                            angle_col=45, annotation_col=group, cluster_rows=TRUE,
  #                            annotation_colors=ann_color, cellwidth=14, cellheight=14)
  #    hmEnrII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
  #                             main=paste0("Genes from ", i, " Enriched Pathway: ", pathway),
  #                             angle_col=45, annotation_col=group, cluster_rows=TRUE,
  #                             annotation_colors=ann_color, cellwidth=14, cellheight=14, )
  #  } else {
  #    hmEnrI[[i]] <- pheatmap(mat=gene.esI[genes, ],
  #                            main=paste0("Genes from ", i, " Enriched Pathway: ", pathway),
  #                            angle_col=45, annotation_col=group, cluster_rows=FALSE,
  #                            annotation_colors=ann_color, cellwidth=14, cellheight=14)
  #    hmEnrII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
  #                             main=paste0("Genes from ", i, " Enriched Pathway: ", pathway),
  #                             angle_col=45, annotation_col=group, cluster_rows=FALSE,
  #                             annotation_colors=ann_color, cellwidth=14, cellheight=14, )
  #  }
    #rowOrder <- hmEnrI[[i]]$tree_row$order
    #hmEnrII[[i]] <- pheatmap(mat=gene.esI[genes, ][rowOrder, ], scale="row",
    #                         main=paste0("Genes from ", i, " Enriched Pathway: ", pathway),
    #                         angle_col=45, annotation_col=group, cluster_rows=FALSE,
    #                         cellwidth=14, cellheight=16, )
  #}

  # plot heatmap for the "core" genes from the intersection of top K enriched pathways and "CORALINE Signature"
  hmCoreEnrI <- vector(mode = "list")
  hmCoreEnrII <- vector(mode = "list")
  for (i in 1:K) {
    pathway <- gseaRes[order(abs(gseaRes$NES), decreasing=TRUE), ]$pathway[i]
    genes <- gseaRes[order(abs(gseaRes$NES), decreasing=TRUE), ]$leadingEdgeGeneNames[[i]]
    genes <- intersect(genes, genes_coraline)
    if (length(genes) > 1) {
      hmCoreEnrI[[i]] <- pheatmap(mat=gene.esI[genes, ],
                                  main=paste0("Core Genes from ", i, " Enriched Pathway: ", pathway),
                                  angle_col=45, annotation_col=group, cluster_rows=TRUE,
                                  annotation_colors=ann_color, cellwidth=14, cellheight=14)
      hmCoreEnrII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
                                   main=paste0("Core Genes from ", i, " Enriched Pathway: ", pathway),
                                   angle_col=45, annotation_col=group, cluster_rows=TRUE,
                                   annotation_colors=ann_color, cellwidth=14, cellheight=14, )
    } else {
      hmCoreEnrI[[i]] <- pheatmap(mat=gene.esI[genes, ],
                                  main=paste0("Core Genes from ", i, " Enriched Pathway: ", pathway),
                                  angle_col=45, annotation_col=group, cluster_rows=FALSE,
                                  annotation_colors=ann_color, cellwidth=14, cellheight=14)
      hmCoreEnrII[[i]] <- pheatmap(mat=gene.esI[genes, ], scale="row",
                                   main=paste0("Core Genes from ", i, " Enriched Pathway: ", pathway),
                                   angle_col=45, annotation_col=group, cluster_rows=FALSE,
                                   annotation_colors=ann_color, cellwidth=14, cellheight=14, )
    }
    #rowOrder <- hmCoreEnrI[[i]]$tree_row$order
    #hmCoreEnrII[[i]] <- pheatmap(mat=gene.esI[genes, ][rowOrder, ], scale="row",
    #                    main=paste0("Core Genes from ", i, " Enriched Pathway: ", pathway),
    #                    angle_col=45, annotation_col=group, cluster_rows=FALSE,
    #                    cellwidth=14, cellheight=16, )
  }

  # save figures
  nes <- arrangeGrob(nes, ncol = 1, nrow = 1)
  ggsave(nes, file=paste0("ESnormalized_", title, "_", Date, ".pdf"), height=80, width=45, limitsize=FALSE)
  #sig <- vector(mode = "list")
  sigCore <- vector(mode = "list")
  #enr <- vector(mode = "list")
  enrCore <- vector(mode = "list")
  for (i in 1:length(hmSigI)) {
  #  sig[[i]] <- arrangeGrob(hmSigI[[i]]$gtable, hmSigII[[i]]$gtable, ncol = 2, nrow = 1)
    sigCore[[i]] <- arrangeGrob(hmCoreSigI[[i]]$gtable, hmCoreSigII[[i]]$gtable, ncol = 2, nrow = 1)
  #  ggsave(sig[[i]], file=paste0("SignificantPathway_", i, "_", type, "_", title, "_", Date, ".pdf"),
  #         height=240, width=50, limitsize=FALSE)
    ggsave(sigCore[[i]], file=paste0("CoreSignificantPathway_", i, "_", type, "_", title, "_", Date, ".pdf"),
           height=H1, width=W1, limitsize=FALSE)
  }
  for (i in 1:length(hmEnrI)) {
  #  enr[[i]] <- arrangeGrob(hmEnrI[[i]]$gtable, hmEnrII[[i]]$gtable, ncol = 2, nrow = 1)
    enrCore[[i]] <- arrangeGrob(hmCoreEnrI[[i]]$gtable, hmCoreEnrII[[i]]$gtable, ncol = 2, nrow = 1)
  #  ggsave(enr[[i]], file=paste0("EnrichedPathway_", i, "_", type, "_", title, "_", Date, ".pdf"),
  #         height=40, width=50, limitsize=FALSE)
    ggsave(enrCore[[i]], file=paste0("CoreEnrichedPathway_", i, "_", type, "_", title, "_", Date, ".pdf"),
           height=H2, width=W2, limitsize=FALSE)
  }

  # return results
  return(list("nes" = nes, "hmCoreSigI" = hmCoreSigI, "hmCoreSigII" = hmCoreSigII,
              "hmCoreEnrI" = hmCoreEnrI, "hmCoreEnrII" = hmCoreEnrII))
}
