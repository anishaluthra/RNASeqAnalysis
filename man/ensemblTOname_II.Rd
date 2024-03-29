% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensemblTOname_II.R
\name{ensemblTOname_II}
\alias{ensemblTOname_II}
\title{Function "ensemblTOname_II"}
\usage{
ensemblTOname_II(geneExp, geneSymbol)
}
\arguments{
\item{geneExp}{a data object with ensembl IDs as rownames}

\item{geneSymbol}{a data frame with all gene annotation information}
}
\value{
a list consisting of a new counts matrix with gene names as rownames.
}
\description{
This function can be called without DESeq, it imports the gene counts matrix based on ensembl IDs and
exports the corresponding matrix based on gene names
}
