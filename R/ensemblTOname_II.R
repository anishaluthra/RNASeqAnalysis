#' Function "ensemblTOname_II"
#'
#' This function can be called without DESeq, it imports the gene counts matrix based on ensembl IDs and
#' exports the corresponding matrix based on gene names
#' @importFrom plyr join
#' @param geneExp a data object with ensembl IDs as rownames
#' @param geneSymbol a data frame with all gene annotation information
#' @return a list consisting of the character vector of gene names and a new counts matrix with gene names as rownames.
#' @export

ensemblTOname_II <- function(geneExp, geneSymbol) {
  ensembl <- data.frame(ensembl_gene_id = row.names(geneExp))
  name <- join(ensembl, geneSymbol[, c("external_gene_name", "ensembl_gene_id")], by="ensembl_gene_id")
  name <- name[!duplicated(name$ensembl_gene_id), ]
  names <- as.character(name$external_gene_name)
  row.names(geneExp) <- names

  # return result
  return(list("name" = name, "geneExp" = geneExp))
}
