#' @title Compute hypergeometric test P value of gene overlap
#' @description Compute hypergeometric test P value of overlap between MAIC genes (or a subset) and a user defined
#'     list of genes If P < user-defined threshold, then a significant overlap can be assumed.
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @param overlap_genes Data frame where the first column lists gene symbols
#' @param n_data_genes Number of MAIC genes to include -- Default = "All" -- chr
#' @param universe Size of gene universe. Either a user defined integer or "HGNC" = All protein
#'     coding genes by HGNC, "FANTOM_L" = genes expressed in lung in FANTOM5, or "FANTOM_S" =
#'     genes expressed in spleen in FANTOM5 -- Default: "HGNC" -- chr
#' @return Hypergeometric test P value
#' @details
#' Input columns for `data_genes` should be (this is the standard output of the MAIC algorithm):
#' * `gene` - HGNC gene symbol - chr
#' * 1...
#' * `uID` - Study unique identifier. Column contains study specific gene score - dbl
#' * n...
#' * `maic_score` - MAIC score for gene - dbl
#' * `contributors` - Studies contributing to MAIC score by method - chr
#'
#' * "HGNC" = All protein coding genes in HGNC = 19220 genes
#' * "FANTOM_L" = Genes expressed in lung in FANTOM5 = 12606 genes
#' * "FANTOM_S" = Genes expressed in spleen in FANTOM5 = 11512 genes
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  overlap_test(ARDSMAICR::data_genes, ARDSMAICR::data_biolitlmine,
#'    n_data_genes = 300, universe = "FANTOM_L")
#'  }
#' }
#' @export
#' @rdname overlap_test
#'
#' @import dplyr
#' @importFrom stats phyper


overlap_test <- function(data_genes, overlap_genes, n_data_genes = "All", universe = 19220) {
  ## Set behaviour for number of maic genes to include

  if (n_data_genes == "All") {
    n_data_genes <- nrow(data_genes)
  } else {
    n_data_genes <- as.numeric(n_data_genes)
  }

  n_condition <- n_data_genes > nrow(data_genes)

  if (n_condition == TRUE) {
    stop("More genes than in MAIC")
  }

  ## Set behaviour for gene universe to compute against

  if (universe == "HGNC") {
    universe <- 19220
  } else if (universe == "FANTOM_L") {
    universe <- 12606
  } else if (universe == "FANTOM_S") {
    universe <- 11512
  } else {
    n_data_genes <- n_data_genes
  }

  ## Subset maic genes based on n_data_genes

  data_genes_slice <- data_genes |>
    dplyr::slice(1:n_data_genes) |>
    dplyr::select(.data$gene)

  ## Extract gene lists

  data_genes_slice <- data_genes_slice |>
    dplyr::pull(.data$gene)

  data_overlap <- overlap_genes |>
    dplyr::select(, 1) %>%
    dplyr::pull()

  ## Find overlap

  overlap <- intersect(data_genes_slice, data_overlap)

  len_overlap <- length(overlap)

  ## hypergeometric test

  len_data_genes <- length(data_genes_slice)

  len_data_overlap <- length(data_overlap)

  hyperp_p <- stats::phyper(len_overlap-1, len_data_overlap, universe-len_data_overlap, len_data_genes, lower.tail = FALSE)

  return(hyperp_p)
}
