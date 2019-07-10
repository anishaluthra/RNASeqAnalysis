#' Function "ensemblTOentrez"
#'
#' This function can be called without DESeq, input a gene expression matrix based on ensembl IDs
#' and output the corresponding matrix based on entrez IDs
#' @importFrom plyr join
#' @param geneExp a data object with ensembl IDs as rownames
#' @param geneSymbol a data frame with all gene annotation information
#' @return a list consisting of the corresponding entrezgene IDs and a new data object with entrezgene ID as rownames
#' @export

ensemblTOentrez <- function(geneExp, geneSymbol) {
  ensembl <- data.frame(ensembl_gene_id = rownames(geneExp))
  entrez <- join(ensembl, geneSymbol[, c("entrezgene_id", "ensembl_gene_id")], by="ensembl_gene_id")
  entrez <- entrez[!duplicated(entrez$ensembl_gene_id), ]
  entrezID <- as.character(entrez$entrezgene)
  rownames(geneExp) <- entrezID

  # return result
  return(list("entrez" = entrez, "geneExp" = geneExp))
}
