#' Function "deseq2"
#'
#' This function Perform Differential Gene Expression Analysis using DESeq2.
#' @import DESeq2
#' @import BiocParallel
#' @param condition a character vector composed of all different group conditions
#' @param n1 a numeric variable indicating the number of samples for one of the group condition
#' @param n2 a numeric variable indicating the number of samples for the other group condition
#' @param N a numeric variable indicating the threshold for filtering low counts genes
#' @param countsmatrix gene expression matrix with raw counts
#' @param M a numeric variable indicating the number of cores for parallel computing
#' @return a list consisting of the DESeq data object and a data frame with all differential gene expression analysis results.
#' @export
#'
deseq2 <- function(condition, n1, n2, N, countsmatrix, M) {
  # set seed
  SEED <- 20190606
  set.seed(SEED)

  # create a data frame coldata to store the group information about all samples
  coldata <- data.frame(condition = c(rep(condition[1], n1),
                                      rep(condition[2], n2)),
                        row.names = colnames(countsmatrix))

  # construct a DESeqDataSet from the countsmatrix and sample information "coldata"
  dds <- DESeqDataSetFromMatrix(countData = countsmatrix,
                                colData = coldata,
                                design = ~ condition)

  # pre-filter the low count genes with reads less than N
  dds_filtered <- dds[rowSums(counts(dds)) > N]

  # set up the factor levels
  dds_filtered$condition <- factor(dds_filtered$condition,
                                   levels = c(condition[1], condition[2]))

  # differencial gene expression analysis
  dds_filtered <- DESeq(dds_filtered, parallel = TRUE, BPPARAM = MulticoreParam(M))

  # generate results
  # default FDR<0.1
  res <- results(dds_filtered, contrast=c("condition", condition[1], condition[2]))

  # return result
  return(list("dds" = dds_filtered, "res" = res))
}
