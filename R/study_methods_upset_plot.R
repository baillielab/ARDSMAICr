#' @title Methods upset plot
#' @description Creates an upset plot of genes, methods, MAIC scores, and their intersections.
#' @param data Data frame in the format of `data_genes`
#' @return An upset plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  methods_upset(data_genes)
#'  }
#' }
#' @export
#' @rdname methods_upset
#'
#' @import dplyr
#' @import stringr
#' @import tibble
#' @import ComplexUpset
#' @import ggbeeswarm
#' @import ggplot2

methods_upset <- function(data) {

  ## Will break if any new method types are identified in a future update

  warning("Only GWAS, TRANSCRIPTOMICS, & PROTEOMICS supported as methods...")

  data_plot <- data |>
    dplyr::mutate(support = dplyr::case_when(
      stringr::str_detect(.data$contributors, "PROTEOMICS") & stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & stringr::str_detect(.data$contributors, "GWAS") ~ "ALL",
      stringr::str_detect(.data$contributors, "PROTEOMICS") & stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & !stringr::str_detect(.data$contributors, "GWAS") ~ "PT",
      stringr::str_detect(.data$contributors, "PROTEOMICS") & stringr::str_detect(.data$contributors, "GWAS") & !stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") ~ "PG",
      stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & stringr::str_detect(.data$contributors, "GWAS") & !stringr::str_detect(.data$contributors, "PROTEOMICS") ~ "TG",
      stringr::str_detect(.data$contributors, "PROTEOMICS") & !stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & !stringr::str_detect(.data$contributors, "GWAS") ~ "P",
      stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & !stringr::str_detect(.data$contributors, "PROTEOMICS") & !stringr::str_detect(.data$contributors, "GWAS") ~ "T",
      stringr::str_detect(.data$contributors, "GWAS") & !stringr::str_detect(.data$contributors, "PROTEOMICS") & !stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") ~ "G",
      TRUE ~ as.character(contributors)
    )) |>
    dplyr::mutate(GWAS = dplyr::case_when(
      stringr::str_detect(.data$support, "ALL") | stringr::str_detect(.data$support, "G") | stringr::str_detect(.data$support, "PG") | stringr::str_detect(.data$support, "TG") ~ "TRUE",
      stringr::str_detect(.data$support, "PT") | stringr::str_detect(.data$support, "P") | stringr::str_detect(.data$support, "T") ~ "FALSE"
    )) |>
    dplyr::mutate(Transcriptomics = dplyr::case_when(
      stringr::str_detect(.data$support, "ALL") | stringr::str_detect(.data$support, "T") | stringr::str_detect(.data$support, "TG") | stringr::str_detect(.data$support, "PT") ~ "TRUE",
      stringr::str_detect(.data$support, "PG") | stringr::str_detect(.data$support, "P") | stringr::str_detect(.data$support, "G") ~ "FALSE"
    )) |>
    dplyr::mutate(Proteomics = dplyr::case_when(
      stringr::str_detect(.data$support, "ALL") | stringr::str_detect(.data$support, "P") | stringr::str_detect(.data$support, "PG") | stringr::str_detect(.data$support, "PT") ~ "TRUE",
      stringr::str_detect(.data$support, "TG") | stringr::str_detect(.data$support, "T") | stringr::str_detect(.data$support, "G") ~ "FALSE"
    )) |>
    dplyr::select(c(.data$gene, .data$GWAS, .data$Transcriptomics, .data$Proteomics, .data$maic_score)) |>
    dplyr::mutate(GWAS = as.logical(.data$GWAS)) |>
    dplyr::mutate(Transcriptomics = as.logical(.data$Transcriptomics)) |>
    dplyr::mutate(Proteomics = as.logical(.data$Proteomics)) |>
    tibble::rownames_to_column(var = "rank")

  methods <- colnames(data_plot)[3:5]

  plot <- ComplexUpset::upset(data_plot,
    methods,
    annotations = list(
      "MAIC Score" = ComplexUpset::upset_annotate("maic_score", ggbeeswarm::geom_quasirandom(na.rm = TRUE, alpha = 0.2))
    ),
    height_ratio = 1,
    width_ratio = 0.1,
    set_sizes = (
      ComplexUpset::upset_set_size(position = "right")
      + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    ),
    sort_intersections_by = "degree"
  )

  return(plot)
}
