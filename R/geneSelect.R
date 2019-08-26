#' Function "geneSelect"
#'
#' This function selects most valuable genes based on different methods
#' @importFrom Biobase exprs
#' @importFrom genefilter varFilter
#' @importFrom CancerSubtypes FSbyMAD
#' @param geneES an ExpressionSet object of gene expression matrix with rows of genes and columns of samples
#' @param type a character variable from one of "CV"(coefficient of variation), "IQR"(interquantile range), "MAD"(madian absolute deviation).
#' @param ratio a numeric variable indicating how much ratio of genes will be selected.
#' @return an ExpressionSet object of selected gene expression matrix.
#' @export
#'
geneSelect <- function(geneES, type, ratio) {
  if (type == "CV") {
    # method I: select most valuable genes based on coefficient of variation
    sd <- sapply(X=1:nrow(exprs(geneES)), FUN=function(x) {
      sd(exprs(geneES)[x, ], na.rm = TRUE)
    })
    mean <- sapply(X=1:nrow(exprs(geneES)), FUN=function(x) {
      mean(exprs(geneES)[x, ], na.rm = TRUE)
    })
    cv <- sd/mean
    names(cv) <- rownames(exprs(geneES))
    # choose the top ratio(%) of genes according to their CV
    genes <- names(cv)[cv > quantile(cv, 1-ratio)]
    geneES1 <- geneES[genes, ]
  }
  else if (type == "IQR") {
    # method II: use IQR to filter the genes and preserve ratio(%)
    geneES1 <- varFilter(geneES, var.func=IQR, var.cutoff=1-ratio, filterByQuantile=TRUE)
  }
  else if (type == "MAD") {
    # method III: use MAD to pick raio(%) most variable genes
    geneExp1 <- FSbyMAD(exprs(geneES), cut.type="topk", value=round(ratio*nrow(exprs(geneES))))
    geneES1 <- geneES[row.names(geneExp1), ]
  }

  # return result
  return(geneES1)
}
