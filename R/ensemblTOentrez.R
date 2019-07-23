#' Function "ensemblTOentrez"
#'
#' This function can be called without DESeq, input a gene expression matrix based on ensembl IDs
#' and output the corresponding matrix based on entrez IDs
#' @importFrom plyr join
#' @param geneExp a data object with ensembl IDs as rownames
#' @param geneSymbol a data frame with all gene annotation information
#' @param condition character vector containing group condition
#' @param n1 the number of one condition in "condition"
#' @param n2 the numbe of the second condition in "condition"
#' @return an ExpressionSet object with matrix based on entrez gene IDs
#' @export

ensemblTOentrez <- function(geneExp, geneSymbol, condition, n1, n2) {
  geneExpI <- as.data.frame(geneExp)
  geneExpI$ensembl_gene_id <- row.names(geneExp)
  geneExpI <- join(geneExpI, geneSymbol[, c("entrezgene_id", "ensembl_gene_id")], by="ensembl_gene_id")
  geneExpI <- aggregate(geneExpI[, 1:ncol(geneExp)], by=list(entrez=geneExpI$entrezgene_id), FUN=mean)
  row.names(geneExpI) <- as.character(geneExpI$entrez)
  geneExpI <- as.matrix(geneExpI)[, -1]

  # create phenotype file
  meta.phe <- data.frame(SampleID = colnames(geneExpI),
                         Group = factor(c(rep(condition[1], n1), rep(condition[2], n2))),
                         row.names = colnames(geneExpI))
  meta_phe <- data.frame(labelDescription=c('SampleID', 'Group'),
                         row.names=c('SampleID', 'Group'))
  gene.phe <- new("AnnotatedDataFrame", data=meta.phe, varMetadata=meta_phe)

  # construct ExpressionSet object from gene.exp and gene.phe
  gene.es <- ExpressionSet(assayData=geneExpI,
                           phenoData=gene.phe,
                           annotation="hgu95av2")

  # return result
  return(gene.es)
}
