#' @title Compute measures of gene overlap
#' @description Compute measures of gene overlap between maic genes or a subset and a user defined
#'     list of genes
#' @param data_genes -- Input in the format of `data_study`
#' @param overlap_genes -- Input where the first column lists gene symbols
#' @param n_data_genes -- Number of maic genes to include -- Default: 'All'
#' @param universe -- Size of gene universe. Either a user defined integer or "HGNC" = All protein
#'     coding genes by HGNC, "FANTOM_L" = genes expressed in lung in FANTOM5, or "FANTOM_S" =
#'     genes expressed in spleen in FANTOM5 -- Default: "HGNC"
#' @return OUTPUT_DESCRIPTION
#' @details "HGNC" = 19220 genes, "FANTOM_L" = 12606 genes, "FANTOM_S" = 11512 genes
#' @examples
#' \dontrun{
#' if(interactive()){
#'  test_gene_overlap(data_genes, data_biolitlmine, n_data_genes = 300, universe = "FANTOM_L")
#'  }
#' }
#' @export
#' @rdname test_gene_overlap
#'
#' @import dplyr
#' @import GeneOverlap

test_gene_overlap <- function(data_genes, overlap_genes, n_data_genes = "All", universe = 19220) {

  ## Set behaviour for number of maic genes to include

  if (n_data_genes == "All") {
    n_data_genes <- nrow(data_genes)
  }
  if (is.numeric(n_data_genes)) {
    n_data_genes <- n_data_genes
  } else {
    stop("More genes than in MAIC")
  }

  ## Set behaviour for gene universe to compute against

  if (universe == "HGNC") {
    universe <- 19220
  }
  if (universe == "FANTOM_L") {
    universe <- 12606
  }
  if (universe == "FANTOM_S") {
    universe <- 11512
  } else {
    n_data_genes <- n_data_genes
  }

  ## Subset maic genes based on n_data_genes

  data_genes_slice <- data_genes |>
    dplyr::slice(1:n_data_genes) |>
    dplyr::select(.data$gene)

  ## Extract gene lists

  data_genes_slice <- data_genes_slice$gene

  data_overlap <- overlap_genes |> dplyr::select(, 1) %>% dplyr::pull()

  ## Create gene overlap object

  go_obj <- GeneOverlap::newGeneOverlap(
    data_genes_slice,
    data_overlap,
    genome.size = universe
  )

  ## Compute gene overlap

  go_obj_res <- GeneOverlap::testGeneOverlap(go_obj)

  return(go_obj_res)
}
