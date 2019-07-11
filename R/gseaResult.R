#' Function "gseaResult"
#' This function stores customized results from GSEA
#'
#' @include entrezTOname.R
#' @importFrom data.table fwrite
#' @param geneSymbol a data frame with all gene annotation information
#' @param gseaRes the result data frame from GSEA
#' @param group a string character indicating the group condition
#' @param gs a string character indicating the name of gene set
#' @param Date Date object obtained from Sys.Date
#' @return a list of sorted GSEA result and significant pathways with different FDR controls.
#' @export

gseaResult <- function(geneSymbol, gseaRes, group, gs, Date) {
  # add gene names of leading edge genes to the result
  # no removing duplicated gene names is each leadin edge gene sets
  LEgenes <- entrezTOname(geneSymbol, gseaRes$leadingEdge)
  gseaRes$leadingEdgeGeneNames <- LEgenes
  fwrite(gseaRes[order(gseaRes$pval), ], file=paste0("gseaResult_", group, "_", gs, "_", Date, ".csv"),
         sep=",", sep2=c(""," ",""))

  # select pathways with FDR<=0.1 and FDR<=0.05
  gseaResI <- gseaRes[which(gseaRes$padj<=0.1), ]
  gseaResI <- gseaResI[order(gseaResI$pval), ]
  gseaResII <- gseaRes[which(gseaRes$padj<=0.05), ]
  gseaResII <- gseaResII[order(gseaResII$pval), ]
  fwrite(gseaResI, file=paste0("gseaResult(FDR=0.1)_", group, "_", gs, "_", Date, ".csv"),
         sep=",", sep2=c(""," ",""))
  fwrite(gseaResII, file=paste0("gseaResult(FDR=0.05)_", group, "_", gs, "_", Date, ".csv"),
         sep=",", sep2=c(""," ",""))

  # return results
  return(list("gseaRes" = gseaRes[order(gseaRes$pval), ], "gseaRes0.1" = gseaResI, "gseaRes0.05" = gseaResII))
}
