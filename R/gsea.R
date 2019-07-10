#' Function "gsea"
#'
#' This function performs GSEA analysis.
#' @importFrom stats aggregate
#' @importFrom Biobase ExpressionSet
#' @importFrom fgsea fgsea
#' @importFrom methods new
#' @param res_anno a data frame with the result of DESeq along with gene annotation
#' @param geneEX normalized gene expression matrix based on ensembl IDs
#' @param geneSet a gene set list, consisting of entrezgene IDs
#' @param condition a string character representing the group condition
#' @param N a numeric variable indicating the number of total samples
#' @param n1 a numeric variable indicating the number of samples for one of the group condition
#' @param n2 a numeric variable indicating the number of samples for the other group condition
#' @param Min a numeric variable representing the minimum size of pathway accepted for GSEA
#' @param Max a numeric variable representing the maximum size of pathway accepted for GSEA
#' @return a list of gene rank based on entrezgene ID, result of GSEA and an ExpressionSet object of normalized counts matrix based on entrezgene IDs.
#' @export

gsea <- function(res_anno, geneEX, geneSet, condition, N, n1, n2, Min, Max) {
  # set seed
  SEED <- 20190606
  set.seed(SEED)

  ## Pre-ranking: rank the gene by "-log10(pval)*sign(fold change)"
  # obtained from DESeq2 analysis
  resI <- res_anno[!is.na(res_anno$entrezgene), ]
  gene.anno <- resI[, c("ensembl_gene_id", "entrezgene_id")]
  gene.anno$metric <- -log10(resI$pvalue)*sign(resI$log2FoldChange)
  geneRanks <- aggregate(gene.anno[, 3], by=list(entrez=gene.anno$entrezgene), FUN=mean)
  geneRank <- geneRanks$x
  names(geneRank) <- as.character(geneRanks$entrez)
  geneRank <- sort(geneRank)

  # create gene annotation file
  gene.anno <- resI[, c("ensembl_gene_id", "entrezgene_id")]
  gene.exp <- as.data.frame(geneEX)
  gene.exp$ensembl_gene_id <- rownames(gene.exp)
  gene.exp <- merge(gene.exp, gene.anno, by="ensembl_gene_id", sort=FALSE)[, -1]
  gene.exp <- aggregate(gene.exp[, 1:N],
                        by=list(entrez=gene.exp$entrezgene), FUN=mean)
  row.names(gene.exp) <- as.character(gene.exp$entrez)
  gene.exp <- as.matrix(gene.exp[, -1])

  # create phenotype file
  meta.phe <- data.frame(SampleID = colnames(gene.exp),
                         Group = factor(c(rep(condition[1], n1), rep(condition[2], n2))),
                         row.names = colnames(gene.exp))
  meta_phe <- data.frame(labelDescription=c('SampleID', 'Group'),
                         row.names=c('SampleID', 'Group'))
  gene.phe <- new("AnnotatedDataFrame", data=meta.phe, varMetadata=meta_phe)

  # construct ExpressionSet object from gene.exp and gene.phe
  gene.es <- ExpressionSet(assayData=gene.exp,
                           phenoData=gene.phe,
                           annotation="hgu95av2")

  # gsea analysis
  gene_fgsea <- fgsea(pathways = geneSet, stats = geneRank,
                      minSize=Min, maxSize=Max, nperm=100000)

  # return result
  # geneES represents patient ID and Quantile normalized between array object
  return(list("geneRank" = geneRank, "geneES" = gene.es, "result_fgsea" = gene_fgsea))
}
