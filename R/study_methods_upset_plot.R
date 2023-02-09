#' @title Methods upset plot
#' @description Creates an upset plot of genes, methods, MAIC scores, and their intersections.
#' @param data Data frame in the format of `ARDSMAICR::data_genes`
#' @return An upset plot
#' @details
#' Input columns for `data_study` should be:
#' * `id` - Integer 1 to n studies - dbl
#' * `First_author` - First author family name - chr
#' * `Article_title` - Article title - chr
#' * `Year` - Year of publication - dbl
#' * `Journal` - Journal - chr
#' * `DOI` - Digital object identifier - dbl
#' * `PMID` - PubMed ID - dbl
#' * `uID` - Unique ID. Format is `First_Author Year PMID` - chr
#' * `Method` - Study method e.g., "GWAS" - chr
#' * `Technology` - Technology used e.g., "Microarray" - chr
#' * `Tissue` - Tissue type sampled e.g., "BALF" - chr
#' * `Cell` - Cell type sampled e.g., "Neutrophils" - chr
#' * `Focus` Study focus e.g., "Susceptibility" - chr
#' * `ARDS_pts` - Total number of patients with ARDS included in study - dbl
#' * `ARDS_definition` - Definition of ARDS used in study - chr
#' * `List_available` - Was the gene list associated with the study retrievable - lgl
#'
#' Input columns for `data_genes` should be (this is the standard output of the MAIC algorithm):
#' * `gene` - HGNC gene symbol - chr
#' * 1...
#' * `uID` - Study unique identifier. Column contains study specific gene score - dbl
#' * n...
#' * `maic_score` - MAIC score for gene - dbl
#' * `contributors` - Studies contributing to MAIC score by method - chr
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  methods_upset(ARDSMAICR::data_genes)
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





#' @title Categories upset plot
#' @description Creates an upset plot of genes, categories, MAIC scores, and their intersections.
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @return An upset plot
#' @details
#' Input columns for `data_study` should be:
#' * `id` - Integer 1 to n studies - dbl
#' * `First_author` - First author family name - chr
#' * `Article_title` - Article title - chr
#' * `Year` - Year of publication - dbl
#' * `Journal` - Journal - chr
#' * `DOI` - Digital object identifier - dbl
#' * `PMID` - PubMed ID - dbl
#' * `uID` - Unique ID. Format is `First_Author Year PMID` - chr
#' * `Method` - Study method e.g., "GWAS" - chr
#' * `Technology` - Technology used e.g., "Microarray" - chr
#' * `Tissue` - Tissue type sampled e.g., "BALF" - chr
#' * `Cell` - Cell type sampled e.g., "Neutrophils" - chr
#' * `Focus` Study focus e.g., "Susceptibility" - chr
#' * `ARDS_pts` - Total number of patients with ARDS included in study - dbl
#' * `ARDS_definition` - Definition of ARDS used in study - chr
#' * `List_available` - Was the gene list associated with the study retrievable - lgl
#'
#' Input columns for `data_genes` should be (this is the standard output of the MAIC algorithm):
#' * `gene` - HGNC gene symbol - chr
#' * 1...
#' * `uID` - Study unique identifier. Column contains study specific gene score - dbl
#' * n...
#' * `maic_score` - MAIC score for gene - dbl
#' * `contributors` - Studies contributing to MAIC score by method - chr
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  categories_upset(ARDSMAICR::data_study, ARDSMAICR::data_genes)
#'  }
#' }
#' @rdname categories_upset
#' @export
#'
#' @importFrom dplyr left_join filter select mutate across group_by summarise ungroup summarize case_when
#' @importFrom tidyselect where
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_glue str_count str_detect
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom ComplexUpset upset upset_annotate upset_set_size
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom ggplot2 theme element_text

categories_upset <- function(data_study, data_genes) {

  data_meta <- dplyr::left_join(data_study, categories, by = "uID")

  categories_df <- data_meta |>
    dplyr::filter(.data$List_available == TRUE) |>
    dplyr::select(c(.data$uID, .data$Method, .data$Category))

  data_categories_count <- data_genes |>
    dplyr::select(-c(.data$maic_score, .data$contributors)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
    tidyr::pivot_longer(
      cols = tidyselect::where(is.logical),
      names_to = "uID",
      values_to = "identified"
    )

  n_categories_gene <- dplyr::left_join(data_categories_count, categories_df, by = "uID", copy = TRUE)

  support_long <- n_categories_gene |>
    dplyr::group_by(.data$gene, .data$Category, .data$Method) |>
    dplyr::summarise(n_lists = sum(.data$identified)) |>
    dplyr::ungroup()

  support <- support_long |>
    dplyr::filter(.data$n_lists >= 1) |>
    dplyr::mutate(Categories = stringr::str_glue("{Method}:{Category}")) |>
    dplyr::select(-c(.data$Category, .data$Method)) |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarize(support = paste(.data$Categories, collapse = ",\ ")) |>
    dplyr::ungroup() |>
    dplyr::mutate(n_categories = 1 + stringr::str_count(.data$support, ","))

  base_df <- dplyr::left_join(data_genes, support, by = "gene")

  df <- base_df |>
    dplyr::mutate(GWAS_genotyping = dplyr::case_when(stringr::str_detect(.data$support, "GWAS:genotyping") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::mutate(GWAS_wes = dplyr::case_when(stringr::str_detect(.data$support, "GWAS:wes") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::mutate(Transcriptomics_microarray = dplyr::case_when(
      stringr::str_detect(.data$support, "Transcriptomics:microarray") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::mutate(Transcriptomics_rnaseq = dplyr::case_when(
      stringr::str_detect(.data$support, "Transcriptomics:RNAseq") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::mutate(Transcriptomics_scrnaseq = dplyr::case_when(
      stringr::str_detect(.data$support, "Transcriptomics:scRNAseq") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::mutate(Proteomics_massspec = dplyr::case_when(
      stringr::str_detect(.data$support, "Proteomics:massspec") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::mutate(Proteomics_other = dplyr::case_when(
      stringr::str_detect(.data$support, "Proteomics:other") ~ TRUE,
      .default = FALSE
    )) |>
    dplyr::select(c(
      .data$gene, .data$GWAS_genotyping, .data$GWAS_wes, .data$Transcriptomics_microarray,
      .data$Transcriptomics_rnaseq, .data$Transcriptomics_scrnaseq, .data$Proteomics_massspec, .data$Proteomics_other,
      .data$maic_score
    )) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "rank")

  category_vars <- colnames(df)[3:9]

  plot <- ComplexUpset::upset(df,
    category_vars,
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
