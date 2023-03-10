% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ards_literature_overlap.R
\name{overlap_count}
\alias{overlap_count}
\title{Count MAIC and alternative gene list overlap}
\usage{
overlap_count(data_genes, data_alternative)
}
\arguments{
\item{data_genes}{Data frame in the format of \code{ARDSMAICR::data_genes}}

\item{data_alternative}{Data frame where the first column is a gene list}
}
\value{
A list
}
\description{
Returns the total number of genes identified in the alternative list, the number of
these genes found in MAIC, and the percentage of alternative list genes found in MAIC.
}
\details{
Input columns for \code{data_genes} should be (this is the standard output of the MAIC algorithm):
\itemize{
\item \code{gene} - HGNC gene symbol - chr
\item 1...
\item \code{uID} - Study unique identifier. Column contains study specific gene score - dbl
\item n...
\item \code{maic_score} - MAIC score for gene - dbl
\item \code{contributors} - Studies contributing to MAIC score by method - chr
}
}
\examples{
\dontrun{
if(interactive()){
 overlap_count(ARDSMAICR::data_genes, ARDSMAICR::data_biolitmine)
 }
}
}
