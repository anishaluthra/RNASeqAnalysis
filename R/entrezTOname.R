#' Function "entrezTOname"
#'
#' This function transfers the entrez gene IDs in geneSet to gene names
#' make sure the input geneSet is a list of pathways with character type entrez gene IDs
#' @importFrom data.table fwrite
#' @importFrom plyr join
#' @param geneSymbol a data frame with all gene annotation information
#' @param geneSet a list of pathways with entrezgene IDs, gene set
#' @return a list of pathways with gene names
#' @export

# Transfer the entrez gene IDs in geneSet to gene names
# make sure the input geneSet is a list of pathways with character type entrez gene IDs

entrezTOname <- function(geneSymbol, geneSet) {
  geneSymbol$entrezgene_id <- as.character(geneSymbol$entrezgene_id)
  geneset <- vector(mode = "list")
  for (i in 1:length(geneSet)) {
    genes <- data.frame(entrezgene_id = geneSet[[i]])
    genes <- join(genes, geneSymbol[, c("entrezgene_id", "external_gene_name")], by="entrezgene_id")
    geneset[[i]] <- as.character(unique(genes$external_gene_name))
  }
  names(geneset) <- names(geneSet)

  # return result
  return("geneset" = geneset)
}
