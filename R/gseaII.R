#' Function "gseaII"
#'
#' This function performs GSEA analysis.
#' @importFrom stats aggregate
#' @importFrom fgsea fgsea
#' @param res_anno a data frame with the result of DESeq along with gene annotation
#' @param geneSet a gene set list, consisting of entrezgene IDs
#' @param attr1 a character variable indicating the based on gene symbols in the input matrix
#' @param attr2 the name of column representing p values in the input matrix
#' @param attr3 the name of column representing log2FoldChange in the input matrix
#' @param Min a numeric variable representing the minimum size of pathway accepted for GSEA
#' @param Max a numeric variable representing the maximum size of pathway accepted for GSEA
#' @return a list of gene rank based on entrezgene ID, result of GSEA.
#' @export

gseaII <- function(res_anno, geneSet, attr1, attr2, attr3, Min, Max) {
  # set seed
  SEED <- 20190606
  set.seed(SEED)

  ## Pre-ranking: rank the gene by "-log10(pval)*sign(fold change)"
  # obtained from DESeq2 analysis
  resI <- res_anno[!is.na(res_anno$entrezgene_id), ]
  gene.anno <- resI[, c(attr1, "entrezgene_id")]
  gene.anno$metric <- -log10(resI[, attr2])*sign(resI[, attr3])
  geneRanks <- aggregate(gene.anno[, 3], by=list(entrez=gene.anno$entrezgene_id), FUN=mean)
  geneRank <- geneRanks$x
  names(geneRank) <- as.character(geneRanks$entrez)
  geneRank <- sort(geneRank)

  # gsea analysis
  gene_fgsea <- fgsea(pathways = geneSet, stats = geneRank,
                      minSize=Min, maxSize=Max, nperm=100000)

  # return result
  # geneES represents patient ID and Quantile normalized between array object
  return(list("geneRank" = geneRank, "result_fgsea" = gene_fgsea))
}
