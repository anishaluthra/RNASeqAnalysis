#' Function "annotation"
#'
#' This function creates gene annotation file using provided file and filter
#' @importFrom biomaRt useMart getBM
#' @param genes a vector of gene IDs such as ensembl ID, gene name
#' @param filter a string that is used as the base for extraction of gene information
#' @param dataSet the name for the dataset to be used
#' @return a data frame with corresponding attributes of the input gene IDs
#' @export

annotation <- function(genes, filter, dataSet = "hsapiens_gene_ensembl") {
  # create gene annotation(provide corresponding gene names and entrez IDs)
  # select a database and a dataset to use
  mymart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = dataSet,
                    host = 'www.ensembl.org')
  # select a dataset and update the Mart object
  #mydataset <- useDataset("hsapiens_gene_ensembl", mart = mymart)

  # construct a main query function to retrieve gene names and entrez gene IDs
  # prepare the ensembl gene IDs we have
  gene_symbol <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id",
                                      "external_gene_name", "description"),
                       filters = filter,
                       values = genes,
                       mart = mymart)

  # return result
  return(gene_symbol)
}
