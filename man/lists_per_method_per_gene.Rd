% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_methods_and_gene_support.R
\name{lists_per_method_per_gene}
\alias{lists_per_method_per_gene}
\title{Count number of gene lists by method by gene}
\usage{
lists_per_method_per_gene(data_study, data_genes)
}
\arguments{
\item{data_study}{Data frame in the format of \code{ARDSMAICR::data_study}}

\item{data_genes}{Data frame in the format of \code{ARDSMAICR::data_genes}}
}
\value{
A tibble
}
\description{
Returns the number of lists each unique gene or SNP is found in stratified by
method (the most concise categorisation of lists).
}
\details{
Input columns for \code{data_study} should be:
\itemize{
\item \code{id} - Integer 1 to n studies - dbl
\item \code{First_author} - First author family name - chr
\item \code{Article_title} - Article title - chr
\item \code{Year} - Year of publication - dbl
\item \code{Journal} - Journal - chr
\item \code{DOI} - Digital object identifier - dbl
\item \code{PMID} - PubMed ID - dbl
\item \code{uID} - Unique ID. Format is \verb{First_Author Year PMID} - chr
\item \code{Method} - Study method e.g., "GWAS" - chr
\item \code{Technology} - Technology used e.g., "Microarray" - chr
\item \code{Tissue} - Tissue type sampled e.g., "BALF" - chr
\item \code{Cell} - Cell type sampled e.g., "Neutrophils" - chr
\item \code{Focus} Study focus e.g., "Susceptibility" - chr
\item \code{ARDS_pts} - Total number of patients with ARDS included in study - dbl
\item \code{ARDS_definition} - Definition of ARDS used in study - chr
\item \code{List_available} - Was the gene list associated with the study retrievable - lgl
}

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
   lists_per_method_per_gene(ARDSMAICR::data_study, ARDSMAICR::data_genes) |>
     dplyr::filter(.data$gene == "IL6")
 }
}
}
