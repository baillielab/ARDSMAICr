#' @title Count gene lists
#' @description Return the total number of included gene lists
#' @param data -- Input in the format of `data_study`
#' @return An integer
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  count_lists(data_study)
#'  }
#' }
#' @export
#' @rdname count_lists

count_lists <- function(data) {
  lists_unique <- unique(data$uID)
  length(lists_unique)
}

#' @title Count gene lists per method type
#' @description Return the total number of gene lists belonging to each method identified
#' @param data -- Input in the format of `data_study`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  lists_per_method(data_study)
#'  }
#' }
#' @export
#' @rdname lists_per_method
#'
#' @import dplyr
#' @importFrom rlang .data

lists_per_method <- function(data) {
  data |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(count = dplyr::n()) |>
    dplyr::mutate(percentage = .data$count / sum(.data$count))
}

#' @title Count genes
#' @description Return the total number of unique genes or SNPs across all gene lists
#' @param data -- Input in the format of `data_genes`
#' @return An integer
#' @examples
#' \dontrun{
#' if(interactive()){
#'  count_genes(data_genes)
#'  }
#' }
#' @export
#' @rdname count_genes

count_genes <- function(data) {
  genes_unique <- unique(data$gene)
  length(genes_unique)
}

#' @title Count genes per gene list
#' @description Return the number of unique genes or SNPs identified in each gene list
#' @param data -- Input in the format of `data_genes`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  print(genes_per_list(data_genes), n = 50)
#'  }
#' }
#' @export
#' @rdname genes_per_list
#'
#' @import dplyr
#' @import tidyselect
#' @import tidyr

genes_per_list <- function(data) {
  data |>
    dplyr::select(-c(.data$maic_score, .data$contributors)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
    dplyr::summarise(dplyr::across(tidyselect::where(is.logical), sum)) |>
    tidyr::pivot_longer(tidyselect::everything(),
      names_to = "list",
      values_to = "n_genes"
    ) |>
    dplyr::arrange(.data$list)
}

#' @title Count genes per method
#' @description Return the number of unique genes or SNPs identified by each method type
#' @param data_study -- Input in the format of `data_study`
#' @param data_genes -- Input in the format of `data_genes`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  genes_per_method(data_study, data_genes)
#'  }
#' }
#' @export
#' @rdname genes_per_method
#'
#' @import dplyr
#' @import tidyselect
#' @import tidyr

genes_per_method <- function(data_study, data_genes) {
  methods <- data_study |>
    dplyr::select(c(.data$uID, .data$Method)) |>
    dplyr::rename(list = .data$uID)

  data_gene_count <- data_genes |>
    dplyr::select(-c(.data$maic_score, .data$contributors)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
    dplyr::summarise(dplyr::across(tidyselect::where(is.logical), sum)) |>
    tidyr::pivot_longer(tidyselect::everything(),
      names_to = "list",
      values_to = "n_genes"
    ) |>
    dplyr::arrange(.data$list)

  n_gene_methods <- dplyr::left_join(data_gene_count, methods, by = "list")

  n_gene_methods |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(n_genes = sum(.data$n_genes)) |>
    dplyr::ungroup()
}

#' @title Count the number of lists each gene or SNP is found in
#' @description For each unique gene or SNP return the number of lists that gene or SNP is
#'     found in
#' @param data -- Input in the format of `data_genes`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  IL6_list_count <- data_genes |>
#'    lists_per_gene() |>
#'    dplyr::filter(.data$gene == "IL6")
#'  print(IL6_list_count)
#'  }
#' }
#' @export
#' @rdname lists_per_gene
#'
#' @import dplyr
#' @import tidyselect

lists_per_gene <- function(data) {
  data |>
    dplyr::select(-c(.data$maic_score, .data$contributors)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
    dplyr::rowwise(.data$gene) |>
    dplyr::summarise(n_lists = sum(dplyr::c_across(tidyselect::where(is.logical)))) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$n_lists))
}

#' @title Count the number of lists per method for each gene or SNP
#' @description Return the number of lists each unique gene or SNP is found in stratified
#'     by method type
#' @param data_study -- Input in the format of `data_study`
#' @param data_genes -- Input in the format of `data_genes`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  IL6_list_per_method_count <- data_genes |>
#'    lists_per_method_per_gene() |>
#'    dplyr::filter(.data$gene == "IL6")
#'  print(IL6_list_per_method_count)
#'  }
#' }
#' @export
#' @rdname lists_per_method_per_gene
#'
#' @import dplyr
#' @import tidyselect
#' @import tidyr

lists_per_method_per_gene <- function(data_study, data_genes) {
  methods <- data_study |>
    dplyr::select(c(.data$uID, .data$Method))

  data_method_count <- data_genes |>
    dplyr::select(-c(.data$maic_score, .data$contributors)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
    tidyr::pivot_longer(cols = tidyselect::where(is.logical),
      names_to = "uID",
      values_to = "identified"
    )

  n_methods_gene <- dplyr::left_join(data_method_count, methods, by = "uID", copy = TRUE)

  n_methods_gene |>
    dplyr::group_by(.data$gene, .data$Method) |>
    dplyr::summarise(n_methods = sum(.data$identified)) |>
    dplyr::ungroup()
}

#' @title Count the number of method types for each gene
#' @description Return the number of method types supporting each unique gene or SNP
#' @param data -- Input in the format of `data_genes`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  IL6_method_support <- data_genes |>
#'    methods_per_gene() |>
#'    dplyr::filter(.data$gene == "IL6")
#'  print(IL6_method_support)
#'  }
#' }
#' @export
#' @rdname methods_per_gene
#'
#' @import dplyr
#' @import stringr
#' @import forcats

methods_per_gene <- function(data) {

  ## Will break if any new method types are identified in a future update

  contributions <- data |>
    dplyr::mutate(support = dplyr::case_when(
      stringr::str_detect(.data$contributors, "PROTEOMICS") & stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & stringr::str_detect(.data$contributors, "GWAS") ~ "ALL",
      stringr::str_detect(.data$contributors, "PROTEOMICS") & stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & !stringr::str_detect(.data$contributors, "GWAS") ~ "PT",
      stringr::str_detect(.data$contributors, "PROTEOMICS") & stringr::str_detect(.data$contributors, "GWAS") & !stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") ~ "PG",
      stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & stringr::str_detect(.data$contributors, "GWAS") & !stringr::str_detect(.data$contributors, "PROTEOMICS") ~ "TG",
      stringr::str_detect(.data$contributors, "PROTEOMICS") & !stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & !stringr::str_detect(.data$contributors, "GWAS") ~ "P",
      stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") & !stringr::str_detect(.data$contributors, "PROTEOMICS") & !stringr::str_detect(.data$contributors, "GWAS") ~ "T",
      stringr::str_detect(.data$contributors, "GWAS") & !stringr::str_detect(.data$contributors, "PROTEOMICS") & !stringr::str_detect(.data$contributors, "TRANSCRIPTOMICS") ~ "G",
      TRUE ~ as.character(.data$contributors)
    )) |>
    dplyr::mutate(support = forcats::as_factor(.data$support)) |>
    dplyr::group_by(.data$gene, .data$support) |>
    dplyr::tally() |>
    dplyr::filter(.data$n > 0) |>
    dplyr::select(c(.data$gene, .data$support)) |>
    dplyr::ungroup() |>
    dplyr::mutate(support = as.character(.data$support))

  contributions |>
    dplyr::mutate(n_methods = dplyr::case_when(
      nchar(.data$support) == 3 ~ 3,
      nchar(.data$support) == 2 ~ 2,
      nchar(.data$support) == 1 ~ 1
    )) |>
    dplyr::arrange(dplyr::desc(.data$n_methods))
}
