#' Function "ensemblTOname_II"
#'
#' This function can be called without DESeq, it imports the gene counts matrix based on ensembl IDs and
#' exports the corresponding matrix based on gene names
#' @importFrom plyr join
#' @param geneExp a data object with ensembl IDs as rownames
#' @param geneSymbol a data frame with all gene annotation information
#' @return a list consisting of a new counts matrix with gene names as rownames.
#' @export

ensemblTOname_II <- function(geneExp, geneSymbol) {
  geneExpI <- as.data.frame(geneExp)
  geneExpI$ensembl_gene_id <- row.names(geneExp)
  geneExpI <- join(geneExpI, geneSymbol[, c("external_gene_name", "ensembl_gene_id")], by="ensembl_gene_id")
  geneExpI <- aggregate(geneExpI[, 1:ncol(geneExp)], by=list(gene_name=geneExpI$external_gene_name), FUN=mean)
  row.names(geneExpI) <- as.character(geneExpI$gene_name)
  geneExpI <- as.matrix(geneExpI)[, -1]

  # return result
  return(list("geneExp" = geneExp))
}
