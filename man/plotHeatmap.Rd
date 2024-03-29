% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHeatmap.R
\name{plotHeatmap}
\alias{plotHeatmap}
\title{Function "plotHeatmap"}
\usage{
plotHeatmap(geneExp, group, myPalette, legendName = "Gene expression",
  method = "ward.D2", distance = "spearman", rowName = FALSE,
  s1 = 15, s2 = 8, s3 = 12, s4 = 6, s5 = 12, s6 = 100, H, W,
  title = "Heatmap", Date)
}
\arguments{
\item{geneExp}{a matrix, rows correspond to genes and columns correspond to samples if it is a gene expression matrix.}

\item{group}{a data frame with each column a different column annotation which needs to be added to the heatmap, rownames correspond to sample ID.}

\item{myPalette}{a named list with color palettes for each annotation, each element in the list should be a named vector.}

\item{legendName}{name of the legend, default is "Gene expression".}

\item{method}{a character type variable indicating the clustering method. Default is "ward.D2".}

\item{distance}{a character type variable indicating the distance used for clustering. Default is "spearman".}

\item{rowName}{a logical variable indicating whether to show the row names. Default is FALSE.}

\item{s1}{a numerical variable indicating the fontsize of legend title. Default is 15.}

\item{s2}{a numerical variable indicating the fontsize of row names. Default is 8.}

\item{s3}{a numerical variable indicating the fontsize of column names. Default is 12.}

\item{s4}{a numerical variable indicating the height and width of grid in defined unit of column legend. Default is 6.}

\item{s5}{a numerical variable indicating the width of the dendrogram for columns and rows. Default is 12.}

\item{s6}{a numerical variable indicating the width of space for rowname and colname annotations. Default is 100.}

\item{H}{a numerical variable indicating the height of saved figure. Default is 15.}

\item{W}{a numerical variable indicating the width of save figure. Default is 15.}

\item{Date}{a Date object obtained from Sys.Date}
}
\value{
a "HeatmapList" class object.
}
\description{
This function plots heat maps using given input matrix and group conditions
}
