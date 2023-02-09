#' @title Table of genes
#' @description Creates a searchable html table of genes ranked by MAIC score.
#' @param data Data frame in the format of `ARDSMAICR::data_genes`
#' @return An html table
#' @details
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
#'  gene_table(ARDSMAICR::data_genes)
#'  }
#' }
#' @export
#' @rdname gene_table
#'
#' @import dplyr
#' @import stringr
#' @import forcats
#' @import tidyselect
#' @import tibble
#' @import htmltools
#' @import reactable
#' @importFrom reactablefmtr espn cell_style pill_buttons icon_assign

gene_table <- function(data) {

  ## Use the contributors column in the standard maic output to identify methods associated with each gene
  ## will break if new methods identified in any subsequent update

  warning("Only GWAS, TRANSCRIPTOMICS, & PROTEOMICS supported as methods...")

  data_info <- data |>
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
    dplyr::mutate(support = forcats::as_factor(.data$support))

  ## Count and tidy the output

  data_info_count <- data_info |>
    dplyr::group_by(.data$gene) |>
    dplyr::mutate(dplyr::across(tidyselect::where(~ is.double(.x)) & !.data$maic_score, ~ as.logical(., na.rm = TRUE))) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$gene) |>
    dplyr::mutate(count = sum(dplyr::across(tidyselect::where(is.logical)))) |>
    dplyr::ungroup() |>
    dplyr::select(c(.data$gene, .data$maic_score, .data$support, .data$count)) |>
    dplyr::mutate(maic_score = round(.data$maic_score, 5)) |>
    tibble::rownames_to_column(var = "rank")

  ## Define colours for each combination of method support

  support_methods_cols <- data_info_count |>
    dplyr::mutate(
      support_methods_cols = dplyr::case_when(
        support == "ALL" ~ "#4daf4a",
        support == "PG" ~ "#118ab2",
        support == "PT" ~ "#e41a1c",
        support == "TG" ~ "#984ea3",
        support == "P" ~ "#ff7f00",
        support == "T" ~ "#a65628",
        TRUE ~ "#ffff33"
      )
    )

  ## Themed html table

  make_table <- function(data) {
    bar_style <- function(width = 1, fill = "#e6e6e6", height = "75%", align = c("left", "right"), color = NULL) {
      align <- match.arg(align)
      if (align == "left") {
        position <- paste0(width * 100, "%")
        image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
      } else {
        position <- paste0(100 - width * 100, "%")
        image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
      }
      list(
        backgroundImage = image,
        backgroundSize = paste("100%", height),
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center",
        color = color
      )
    }
    htmltools::tags$link(
      href = "https://fonts.googleapis.com/css?family=Work+Sans:400,600,700&display=fallback",
      rel = "stylesheet"
    )

    reactable::reactable(
      data,
      theme = reactablefmtr::espn(),
      fullWidth = FALSE,
      style = list(fontFamily = "Work Sans, sans-serif", fontSize = "30px"),
      showSortIcon = TRUE,
      searchable = TRUE,
      language = reactable::reactableLang(
        searchPlaceholder = "SEARCH FOR A GENE..."
      ),
      showSortable = TRUE,
      defaultSorted = list(maic_score = "desc", rank = "asc"),
      showPageSizeOptions = TRUE,
      pageSizeOptions = c(25, 50, 100, 500),
      defaultPageSize = 50,
      highlight = TRUE,
      borderless = FALSE,
      outlined = TRUE,
      resizable = TRUE,
      columns = list(
        rank = reactable::colDef(
          name = "Rank",
          sortable = FALSE,
          minWidth = 75
        ),
        gene = reactable::colDef(
          name = "Gene Name",
          sortable = FALSE,
          minWidth = 175,
          style = reactablefmtr::cell_style(
            horizontal_align = "left",
            font_weight = "900",
            font_size = "14px"
          )
        ),
        maic_score = reactable::colDef(
          name = "MAIC score",
          minWidth = 175,
          style = function(value) {
            bar_style(width = value / max(data$maic_score), fill = "#2c5e77", color = "#fff", align = "right")
          }
        ),
        support = reactable::colDef(
          name = "Supporting methodologies",
          sortable = FALSE,
          minWidth = 300,
          align = "center",
          cell = reactablefmtr::pill_buttons(data, color_ref = "support_methods_cols", opacity = 0.7),
          style = list(borderRight = "1px solid #777")
        ),
        count = reactable::colDef(
          name = "Total number of lists",
          minWidth = 175,
          align = "left",
          cell = reactablefmtr::icon_assign(
            data,
            fill_color = "#99d8c9",
            empty_opacity = 0
          )
        ),
        support_methods_cols = reactable::colDef(
          show = FALSE
        )
      )
    )
  }

  g_table <- make_table(support_methods_cols)

  return(g_table)
}
