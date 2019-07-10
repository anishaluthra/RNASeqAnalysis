#' Function "geneExp"
#'
#' This function constructs expression set for GSEA and ssGSEA
#' @importFrom stats aggregate
#' @importFrom Biobase ExpressionSet
#' @importFrom methods new
#' @param res_anno a data frame with the result of DESeq along with gene annotation
#' @param geneEX normalized gene expression matrix based on ensembl IDs
#' @param condition a string character representing the group condition
#' @param N a numeric variable indicating the number of total samples
#' @param n1 a numeric variable indicating the number of samples for one of the group condition
#' @param n2 a numeric variable indicating the number of samples for the other group condition
#' @return a gene expression object with matrix based on entrezgene IDs
#' @export

geneEXP <- function(res_anno, geneEX, condition, N, n1, n2) {
  # create gene annotation file
  resI <- res_anno[!is.na(res_anno$entrezgene), ]
  gene.anno <- resI[, c("ensembl_gene_id", "entrezgene_id")]
  gene.exp <- as.data.frame(geneEX)
  gene.exp$ensembl_gene_id <- row.names(gene.exp)
  gene.exp <- merge(gene.exp, gene.anno, by="ensembl_gene_id", sort=FALSE)[, -1]
  gene.exp <- aggregate(gene.exp[, 1:N],
                        by=list(entrez=gene.exp$entrezgene_id), FUN=mean)
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

  # return result
  return(list("geneES" = gene.es))
}
