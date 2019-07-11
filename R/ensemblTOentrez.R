#' Function "ensemblTOentrez"
#'
#' This function can be called without DESeq, input a gene expression matrix based on ensembl IDs
#' and output the corresponding matrix based on entrez IDs
#' @importFrom plyr join
#' @param geneExp a data object with ensembl IDs as rownames
#' @param geneSymbol a data frame with all gene annotation information
#' @return a list consisting of a new data object with entrezgene ID as rownames
#' @export

ensemblTOentrez <- function(geneExp, geneSymbol) {
  geneExpI <- as.data.frame(geneExp)
  geneExpI$ensembl_gene_id <- row.names(geneExp)
  geneExpI <- join(geneExpI, geneSymbol[, c("entrezgene_id", "ensembl_gene_id")], by="ensembl_gene_id")
  geneExpI <- aggregate(geneExpI[, 1:ncol(geneExp)], by=list(entrez=geneExpI$entrezgene_id), FUN=mean)
  row.names(geneExpI) <- as.character(geneExpI$entrez)
  geneExpI <- as.matrix(geneExpI)[, -1]

  # return result
  return(list("geneExp" = geneExpI))
}
