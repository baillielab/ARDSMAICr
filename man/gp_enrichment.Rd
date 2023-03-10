% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_enrichment.R
\name{gp_enrichment}
\alias{gp_enrichment}
\title{Functional enrichment using gprofiler}
\usage{
gp_enrichment(
  data,
  p_threshold = 0.01,
  source = c("GO", "KEGG", "REAC", "WP", "HPA")
)
}
\arguments{
\item{data}{Data frame in the format of \code{ARDSMAICR::data_genes}}

\item{p_threshold}{Threshold P value -- Default = 0.01}

\item{source}{Database -- any combination of GO, GO:MF, GO:CC, GO:BP, KEGG, REAC, WP, or HPA
-- Default = c("GO", "KEGG", "REAC", "WP", "HPA")}
}
\value{
Enrichment result and associated metadata
}
\description{
Performs functional enrichment of standard MAIC output against a range of databases.
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

Databases:
\itemize{
\item GO = Gene Ontology (MF = molecular function; CC = cellular component, BP = biological process)
\item KEGG = Kyoto Encyclopaedia of Genes and Genomes
\item REAC = Reactome
\item WP = WikiPathways
\item HPA = Human Protein Atlas
}
}
\examples{
\dontrun{
if(interactive()){
 gp_enrichment(ARDSMAICR::data_genes, p_threshold = 0.001, source = c("REAC", "WP"))
 }
}
}
