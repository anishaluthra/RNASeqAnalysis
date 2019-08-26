#' Function "nameTOentrez"
#'
#' This function transfers the name of genes in geneSet to entrez gene IDs,
#' make sure the input geneSet is a list of pathways with "character" type of gene names
#' @importFrom magrittr %>%
#' @param geneSymbol a data frame with all gene annotation information
#' @param geneSet a list of gene set based on gene names
#' @return a list type gene set based on entrezgene IDs.
#' @export

nameTOentrez <- function(geneSymbol, geneSet) {
  geneSymbol$external_gene_name <- as.character(geneSymbol$external_gene_name)
  geneset <- vector(mode = "list")
  for (i in 1:length(geneSet)) {
    genes <- data.frame(external_gene_name = geneSet[[i]])
    genes$external_gene_name <- as.character(genes$external_gene_name)
    genes <- left_join(genes, geneSymbol[, c("entrezgene_id", "external_gene_name")])
    geneset[[i]] <- as.character(genes$entrezgene_id)
    geneset[[i]] <- geneset[[i]][!is.na(geneset[[i]])] %>% unique()
  }
  names(geneset) <- names(geneSet)

  # remove empty pathway
  geneset <- geneset[lapply(geneset, length) > 0]

  # return result
  return("geneset" = geneset)
}
