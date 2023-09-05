#' Study metadata - ARDS MAIC Update 1.0.0
#'
#' A description of each gene list included in the systematic review and ARDS MAIC
#'
#'
#' @format ## `data_study`
#' A data frame with 44 rows and 16 columns:
#' \describe{
#'   \item{id}{integer}
#'   \item{First_author}{First author's surname}
#'   \item{Article_title}{Title}
#'   \item{Year}{Publication year}
#'   \item{Journal}{Journal}
#'   \item{DOI}{Digital Object Identifier}
#'   \item{PMID}{Pubmed ID}
#'   \item{uID}{Unique identifier}
#'   \item{Method}{Method used to create list}
#'   \item{Technology}{Specific technology used}
#'   \item{Tissue}{Tissue samples derived from}
#'   \item{Cell}{Cell type sampled}
#'   \item{Focus}{Study objective}
#'   \item{ARDS_pts}{Total number of ARDS patients included in study}
#'   \item{ARDS_definition}{ARDS definition used in study}
#'   \item{List_available}{Is a gene list retrievable}
#'   ...
#' }
#' @source Own
"data_study"

#' Overall ARDS MAIC results - ARDS MAIC Update 1.0.0
#'
#' ARDS MAIC output for all gene lists
#'
#'
#' @format ## `data_genes`
#' A data frame with 7088 rows and 46 columns:
#' \describe{
#'   \item{gene}{HGNC symbol}
#'   \item{contributors}{Studies contributing to MAIC score}
#'   ...
#' }
#' @source Own
"data_genes"

#' Overall gene list contributions - ARDS MAIC Update 1.0.0
#'
#' Gene list information content and information contribution
#'
#'
#' @format ## `data_contributions`
#' A data frame with 43 rows and 5 columns:
#' \describe{
#'   \item{study}{Gene list}
#'   \item{raw_information}{Raw information}
#'   \item{information}{Information}
#'   \item{raw_contribution}{Raw contribution}
#'   \item{contribution}{Contribution}
#'   ...
#' }
#' @source Own
"data_contributions"

#' BioLitMine results - ARDS MAIC Update 1.0.0
#'
#' Genes associated with the MeSH heading 'Respiratory Distress syndrome, Acute' obtained on the
#'     4th January, 2023
#'
#'
#' @format ## `data_biolitmine`
#' A data frame with 271 rows and 3 columns:
#' \describe{
#'   \item{Gene}{HGNC symbol}
#'   \item{Last_10_PMIDS}{The PubMed IDs of the last 10 associated publications}
#'   \item{Publication_count}{The total number of associated publications}
#'   ...
#' }
#' @source \url{https://www.flyrnai.org/tools/biolitmine/web/}
"data_biolitmine"

#' SARS-CoV-2 MAIC genes - ARDS MAIC Update 1.0.0
#'
#' Genes included in a MAIC of SARS-CoV-2 genome-wide studies
#'
#'
#' @format ## `data_covidmaic`
#' A data frame with 5418 rows and 1 column:
#' \describe{
#'   \item{gene}{HGNC symbol}
#'   ...
#' }
#' @source \url{https://doi.org/10.1038/s41598-020-79033-3}
"data_covidmaic"
