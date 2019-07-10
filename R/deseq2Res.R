#' Function "deseq2Res"
#'
#' This function saves the differential gene expression analysis results and retrieves raw and normalized counts matrices.
#' @import DESeq2
#' @importFrom data.table fwrite
#' @param countsmatrix raw counts matrix
#' @param res the result data frame from DESeq2
#' @param geneSymbol a data frame with all gene annotation information
#' @param index a character used to merge the result and gene annotation file
#' @param group a string character indicating the group condition
#' @param dds a DESeq data object retrieved from DESeq2 used for accessing normalized counts matrix
#' @param Date a Date object obtained from Sys.Date
#' @return a list consisting of the sorted DESeq result, differentially expressed gene lists, raw and normalized counts matrices.
#' @export

deseq2Res <- function(countsmatrix, res, geneSymbol, index, group, dds, Date) {
  # sort the result by adjusted p values
  res_sort <- res[order(res$padj), ]

  # combine results file and gene annotation file
  res_anno <- as.data.frame(res)
  res_anno$ensembl_gene_id <- rownames(res)
  res_anno <- merge(res_anno, geneSymbol, by=index)

  # pick differentially expressed genes at FDR=0.1 and FDR=0.05
  genesDE <- res_anno[which(res_anno$padj <= 0.1), ]
  genesDE <- genesDE[order(genesDE$padj), ]
  genesDE_I <- res_anno[which(res_anno$padj <= 0.05), ]
  genesDE_I <- genesDE_I[order(genesDE_I$padj), ]
  fwrite(genesDE, file=paste0("DEgenes_", group, "_FDR=0.1_", Date, ".csv"))
  fwrite(genesDE_I, file=paste0("DEgenes_", group, "_FDR=0.05_", Date, ".csv"))

  # remove duplicate ensembls
  genesDE_rd <- genesDE[!duplicated(genesDE$ensembl_gene_id), ]
  genesDE_I_rd <- genesDE_I[!duplicated(genesDE_I$ensembl_gene_id), ]
  fwrite(genesDE_rd, file=paste0("DEgenes(no duplicates)_", group, "_FDR=0.1_", Date, ".csv"))
  fwrite(genesDE_I_rd, file=paste0("DEgenes(no duplicates)_", group, "_FDR=0.05_", Date, ".csv"))

  # sort the result by adjusted p values
  res_anno_ordered <- res_anno[order(res_anno$padj), ]
  fwrite(res_anno_ordered, file=paste0("DESeq2result_", group, "_", Date, ".csv"))

  # access the normalized and unnormalized counts matrix (normalized by median ratio method)
  dds_nor <- estimateSizeFactors(dds)
  dds_nor <- counts(dds_nor, normalized=TRUE)
  dds_unor <- counts(estimateSizeFactors(dds), normalized=FALSE)
  raw <- countsmatrix[which(rownames(countsmatrix) %in% rownames(dds_unor)), ]

  # return results
  return(list("res_sort" = res_sort, "res_anno" = res_anno, "genesDE0.1" = genesDE,
              "genesDE0.05" = genesDE_I, "genesDE0.1_rd" = genesDE_rd, "genesDE0.05_rd" = genesDE_I_rd,
              "res_anno_ordered" = res_anno_ordered, "dds_nor" = dds_nor,
              "dds_unor" = dds_unor, "raw" = raw))
}
