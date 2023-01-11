#' @title Functional enrichment using gprofiler
#' @description Performs functional enrichment of standard MAIC output against a range of databases.
#' @param data Data frame in the format of `data_genes`
#' @param p_threshold Threshold P value -- Default = 0.01
#' @param source Database -- any combination of GO, GO:MF, GO:CC, GO:BP, KEGG, REAC, WP, or HPA
#'     -- Default = c("GO", "KEGG", "REAC", "WP", "HPA")
#' @return Enrichment result and associated metadata
#' @details
#' Input columns for `data_genes` should be (this is the standard output of the MAIC algorithm):
#' * `gene` - HGNC gene symbol - chr
#' * 1...
#' * `uID` - Study unique identifier. Column contains study specific gene score - dbl
#' * n...
#' * `maic_score` - MAIC score for gene - dbl
#' * `contributors` - Studies contributing to MAIC score by method - chr
#'
#' Databases:
#' * GO = Gene Ontology (MF = molecular function; CC = cellular component, BP = biological process)
#' * KEGG = Kyoto Encyclopaedia of Genes and Genomes
#' * REAC = Reactome
#' * WP = WikiPathways
#' * HPA = Human Protein Atlas
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  gp_enrichment(data_genes, p_threshold = 0.001, source = c("REAC", "WP"))
#'  }
#' }
#' @export
#' @rdname gp_enrichment
#'
#' @import dplyr
#' @import gprofiler2

gp_enrichment <- function(data, p_threshold = 0.01, source = c("GO", "KEGG", "REAC", "WP", "HPA")) {

  ## Extract gene names from standard maic output

  gene_names <- data |>
    dplyr::pull(.data$gene)

  ## Conduct functional enrichment for H. sapiens against GO, KEGG, Reactome, WikiPathways, and
  ## the Human Protein Atlas add a P threshold of 0.01

  gostres <- gprofiler2::gost(
    query = gene_names,
    organism = "hsapiens",
    ordered_query = TRUE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = p_threshold,
    correction_method = "g_SCS",
    domain_scope = "annotated",
    custom_bg = NULL,
    numeric_ns = "",
    sources = source,
    as_short_link = FALSE
  )

  return(gostres)
}

#' @title Functional enrichment plot
#' @description Creates a Manhattan plot of functional enrichment results.
#' @param gp_result The return object of `gp_enrichment()`
#' @param is_interactive Interactive plot -- Default = TRUE
#' @return A Manhattan plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  gp_enrichment(data_genes, p_threshold = 0.001, source = c("REAC", "WP")) |>
#'  gp_plot(is_interactive = FALSE)
#'  }
#' }
#' @export
#' @rdname gp_plot
#'
#' @import gprofiler2

gp_plot <- function(gp_result, is_interactive = TRUE) {

  plot <- gprofiler2::gostplot(gp_result, capped = TRUE, interactive = is_interactive)

  return(plot)
}

#' @title Tidy functional enrichment results
#' @description Returns a tidy tibble of the gprofiler output.
#' @param gp_result The return object of `gp_enrichment()`
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  gp_enrichment(data_genes, p_threshold = 0.001, source = c("REAC", "WP")) |>
#'  gp_tidy_results()
#'  }
#' }
#' @export
#' @rdname gp_tidy_results
#'
#' @import dplyr
#' @import tidyselect
#' @import tibble
#' @import stringr
#' @import forcats

gp_tidy_results <- function(gp_result) {

  ## Extract result

  result <- gp_result$result

  ## Tidy

  result %>%
    dplyr::arrange(dplyr::desc(.data$recall)) |>
    dplyr::select(c(.data$p_value, .data$term_size, .data$query_size, .data$intersection_size, .data$precision, .data$recall, .data$term_id, .data$source, .data$term_name)) |>
    dplyr::relocate(c("term_name", "source", "term_id"), .before = "p_value") |>
    dplyr::relocate("intersection_size", .before = "term_size") |>
    dplyr::relocate(c("recall", "precision"), .before = "query_size") |>
    dplyr::relocate("p_value", .after = tidyselect::last_col()) |>
    tibble::rownames_to_column("Rank") |>
    dplyr::mutate(term_name = stringr::str_to_sentence(.data$term_name)) |>
    dplyr::mutate(source = forcats::as_factor(.data$source)) |>
    dplyr::mutate(p_value = formatC(.data$p_value, format = "e", digits = 2))
}

#' @title Functional enrichment results table
#' @description Creates an interactive html table of the tidy output of gprofiler.
#' @param gp_tidy_result The return object of `gp_tidy_results()`
#' @return An html table
#' @examples
#' \dontrun{
#' if(interactive()){
#'  gp_enrichment(data_genes, p_threshold = 0.001, source = c("REAC", "WP")) |>
#'  gp_tidy_results() |>
#'  gp_table()
#'  }
#' }
#' @export
#' @rdname gp_table
#'
#' @import htmltools
#' @import reactable
#' @importFrom reactablefmtr espn

gp_table <- function(gp_tidy_result) {

  ## Set style sheet

  htmltools::tags$link(
    href = "https://fonts.googleapis.com/css?family=Work+Sans:400,600,700&display=fallback",
    rel = "stylesheet"
  )

  ## Html table

  table <- reactable::reactable(
    gp_tidy_result,
    theme = reactablefmtr::espn(),
    fullWidth = TRUE,
    style = list(fontFamily = "Work Sans, sans-serif", fontSize = "30px"),
    showSortIcon = TRUE,
    searchable = TRUE,
    language = reactable::reactableLang(
      searchPlaceholder = "PATHWAY..."
    ),
    showSortable = TRUE,
    defaultSorted = list(Rank = "desc", recall = "asc"),
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(25, 50, 100, 500),
    defaultPageSize = 50,
    highlight = TRUE,
    borderless = FALSE,
    outlined = TRUE,
    resizable = TRUE,
    columns = list(
      Rank = reactable::colDef(show = FALSE),
      term_name = reactable::colDef(
        name = "Pathway",
        cell = function(value, index) {
          term_id <- gp_tidy_result$term_id[index]
          term_id <- if (!is.na(term_id)) term_id else "Unknown"
          htmltools::div(
            htmltools::div(style = list(fontWeight = 600, fontsize = 10), value),
            htmltools::div(style = list(fontsize = 10), term_id),
          )
        },
        minWidth = 300
      ),
      source = reactable::colDef(show = FALSE),
      term_id = reactable::colDef(show = FALSE),
      intersection_size = reactable::colDef(name = "Number of pathway genes in MAIC"),
      term_size = reactable::colDef(name = "Total number of genes in pathway"),
      recall = reactable::colDef(
        name = "Proportion of total pathway genes in MAIC",
        format = reactable::colFormat(digits = 3)
      ),
      precision = reactable::colDef(
        name = "Proportion of MAIC genes in pathway",
        format = reactable::colFormat(digits = 3)
      ),
      query_size = reactable::colDef(
        name = "MAIC rank limit for pathway",
        format = reactable::colFormat(digits = 0)
      ),
      p_value = reactable::colDef(
        name = "P value for pathway",
        align = "right"
      )
    )
  )

  return(table)
}

#' @title Count enriched pathways
#' @description Returns the total number of enriched pathways or the number of enriched pathways per
#'     database.
#' @param gp_result The return object of `gp_enrichment()`
#' @param by_pathway Count by pathway database -- Default = FALSE
#' @return Either an integer or a tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  gp_enrichment(data_genes, p_threshold = 0.001, source = c("REAC", "WP")) |>
#'  gp_count_pathways()
#'  }
#' }
#' @export
#' @rdname gp_count_pathways
#'
#' @import dplyr
#' @import forcats

gp_count_pathways <- function(gp_result, by_pathway = FALSE){

  ## Extract result

  result <- gp_result$result

  ## Set behaviour for by_pathway = TRUE or FALSE

  if (by_pathway) {

    ## TRUE

    result |>
      dplyr::mutate(source = forcats::as_factor(.data$source)) |>
      dplyr::count(.data$source, sort = TRUE)

  } else {

    ## FALSE

    result |>
      dplyr::mutate(term_id = forcats::as_factor(.data$term_id)) |>
      dplyr::group_by(.data$term_id) |>
      dplyr::summarize(distinct_terms = dplyr::n_distinct(.data$term_id)) |>
      nrow()
  }
}






