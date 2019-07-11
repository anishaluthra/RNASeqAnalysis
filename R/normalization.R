#' Function "normalization"
#'
#' This function performs gene expression normalization.
#' @importFrom cqn cqn
#' @param res unnormalized gene expresison matrix based on ensembl IDs
#' @param geneLGC a data frame with gene length and GC content information
#' @return a list consisting of normalized gene expression matrix.
#' @export

normalization <- function(res, geneLGC) {
  # set seed
  SEED <- 20190606
  set.seed(SEED)

  ## Gene expression normalization(use conditional quantile normalization)
  gccontent <- rep(0, nrow(res))
  gccontent <- sapply(X=1:nrow(res), FUN=function(x) {
    gccontent[x] <- geneLGC$`Gene % GC content`[which(geneLGC$ensemble_gene_id==rownames(res)[x])]
  })
  glength <- rep(0, nrow(res))
  glength <- sapply(X=1:nrow(res), FUN=function(x) {
    glength[x] <- geneLGC$`Transcript length (including UTRs and CDS)`[which(geneLGC$ensemble_gene_id==rownames(res)[x])]
  })
  cqn <- cqn(counts=res, x=gccontent, lengths=glength, lengthMethod="smooth")
  geneExp <- cqn$y + cqn$offset

  # return result
  return(list("cqn" = cqn, "geneExp" = geneExp))
}
