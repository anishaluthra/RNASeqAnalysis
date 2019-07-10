#' Function "dataPreprocess"
#'
#' This function preprocesses clinical data and gene expression matrix.
#' @param dataInfo a data frame which contains clinical information of data
#' @param dataMatrix a matrix consisting of gene expression counts
#' @param pcg_ense a character vector of protein-coding ensembl gene IDs
#' @param condition a character vector composed of all different group conditions
#' @param attr1 a string character indicating the column name for "Sample ID" in the clinical data
#' @param attr2 a string character indicating the desired group name in the clinical data
#' @param type1 a string character representing one group condition
#' @param type2 a string character representing the other group condition
#' @param Order a logical variable indicating whether or not to use decreasing order to resort the data information
#' @return a list consisting of numbers of each group condition, processed clinical data information and ready-to-use gene expression matrix.
#' @export

dataPreprocess <- function(dataInfo, dataMatrix, pcg_ense, condition, attr1, attr2, type1, type2, Order) {
  # select the desired columns from dataInfo
  Info <- dataInfo[, c(attr1, attr2)]
  Info <- Info[which(Info[, attr2] %in% c(type1, type2)), ]

  # create a Type column to store the group information of each sample
  Info$Type <- rep("type", dim(Info)[1])

  # figure out which samples belong to type1 and which belong to type2
  indices1 <- which(Info[, attr2] == type1)
  indices2 <- which(Info[, attr2] == type2)
  N1 <- length(indices1)
  N2 <- length(indices2)
  Info$Type[indices1] <- lapply(X=1:N1, FUN=function(x){paste0(condition[1], x)})
  Info$Type[indices2] <- lapply(X=1:N2, FUN=function(x){paste0(condition[2], x)})
  Info_sorted <- Info[order(Info[, attr2], decreasing = Order), ]
  Info_sorted$Group <- c(rep(condition[1], N1), rep(condition[2], N2))

  # construct counts matrix
  countsmatrix <- dataMatrix[, -1]
  rownames(countsmatrix) <- as.character(dataMatrix[, 1])
  class(countsmatrix) <- "integer"

  # remove the version information ".*" from the gene IDs
  rownames(countsmatrix) <- lapply(X=1:nrow(countsmatrix), FUN=function(x){
    strsplit(rownames(countsmatrix)[x], "[.]")[[1]][1]
  })

  # construct counts matrix for only protein-coding related ensembl ids
  countsmatrix <- countsmatrix[which(rownames(countsmatrix) %in% pcg_ense), ]

  # filter out those samples not included in Info_sorted
  countsmatrix <- countsmatrix[, which(colnames(countsmatrix) %in% as.character(Info_sorted[, attr1]))]

  # create a sorted counts matrix by its group conditions
  countsmatrixI <- countsmatrix
  countsmatrixI <- sapply(X=1:ncol(countsmatrixI), FUN=function(x)
  {countsmatrixI[, x] <- countsmatrixI[, which(colnames(countsmatrixI)==as.character(Info_sorted[, attr1][x]))]})
  colnames(countsmatrixI) <- as.character(Info_sorted[, attr1])

  # return results
  return(list("N1" = N1, "N2" = N2, "countsmatrix" = countsmatrixI, "Info" = Info_sorted))
}
