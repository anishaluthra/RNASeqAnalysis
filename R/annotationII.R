#' Function "annotationII"
#'
#' This function transfers genes from one species to another
#' @importFrom biomaRt useMart getLDS
#' @param genes a vector of gene IDs such as ensembl ID, gene name
#' @param filter a string that is used as the base for extraction of gene information
#' @return a data frame with corresponding attributes of the input gene IDs
#' @export

annotationII <- function(genes, dataSet1="hsapiens_gene_ensembl", dataSet2="mmusculus_gene_ensembl",
                         filter1=c("mgi_symbol"), filter2=c("external_gene_name", "ensembl_gene_id", "entrezgene_id")) {
  # select a database and a dataset to use
  mart1 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = dataSet1,
                   host = 'www.ensembl.org')
  mart2 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = dataSet2,
                   host = 'www.ensembl.org')

  # construct a main query function to retrieve gene names and entrez gene IDs
  # prepare the ensembl gene IDs we have
  geneSymbol <- getLDS(attributes = filter1, filters = filter1, values = genes , mart = mart2,
                       attributesL = filter2, martL = mart1, uniqueRows=T)

  # return result
  return(geneSymbol)
}
