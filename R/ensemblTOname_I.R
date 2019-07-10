#' Function "ensemblTOname_I"
#'
#' This function should be called after DESeq, it transfers counts matrix based on ensembl IDs to matrix
#' based on gene names and takes average for duplicated gene names
#' @importFrom stats aggregate
#' @importFrom Biobase ExpressionSet
#' @importFrom methods new
#' @param res_anno a data frame with the result of DESeq along with gene annotation
#' @param raw the raw counts matrix
#' @param N a numeric value representing the number of samples
#' @param condition a string character representing the group condition
#' @param n1 a numeric variable indicating the number of samples for one of the group condition
#' @param n2 a numeric variable indicating the number of samples for the other group condition
#' @return a ExpressionSet object of counts matrix based on gene names with group condition information
#' @export

ensemblTOname_I <- function(res_anno, raw, N, condition, n1, n2) {
  # create gene annotation file
  gene.anno <- res_anno[, c("ensembl_gene_id", "external_gene_name")]
  gene.exp <- as.data.frame(raw)
  gene.exp$ensembl_gene_id <- rownames(gene.exp)
  gene.exp <- merge(gene.exp, gene.anno, by="ensembl_gene_id", sort=FALSE)[, -1]
  gene.exp <- aggregate(gene.exp[, 1:N],
                        by=list(gene_name=gene.exp$external_gene_name), FUN=mean)
  row.names(gene.exp) <- gene.exp$gene_name
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
  return(gene.es)
}
