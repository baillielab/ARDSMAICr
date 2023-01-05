#' @title Create a table of study designs
#' @description Create an interactive html table of study designs
#' @param data -- Input in the format of `data_study`
#' @return An html table
#' @examples
#' \dontrun{
#' if(interactive()){
#'  design_table(data_study)
#'  }
#' }
#' @export
#' @rdname design_table
#'
#' @import dplyr
#' @import reactable
#' @importFrom reactablefmtr espn cell_style pill_buttons icon_assign

design_table <- function(data) {

  ## Subset table

  data <- data |>
    dplyr::select(c(.data$First_author, .data$PMID, .data$Year, .data$Focus, .data$ARDS_definition, .data$ARDS_pts)) |>
    dplyr::distinct(.data$PMID, .keep_all = TRUE) # will need preprint tag to change if > 1 at next update

  ## Themed html table

  table <- reactable::reactable(
    data,
    theme = reactablefmtr::espn(),
    fullWidth = TRUE,
    style = list(fontFamily = "Work Sans, sans-serif", fontSize = "20px"),
    showSortIcon = TRUE,
    showSortable = TRUE,
    searchable = FALSE,
    pagination = FALSE,
    defaultSorted = list(First_author = "asc", Year = "desc"),
    highlight = TRUE,
    borderless = TRUE,
    striped = TRUE,
    outlined = TRUE,
    resizable = TRUE,
    wrap = FALSE,
    columns = list(
      First_author = reactable::colDef(name = "Study", sticky = "left", minWidth = 75),
      PMID = reactable::colDef(name = "PMID", sticky = "left", minWidth = 75),
      Year = reactable::colDef(name = "Year", sticky = "left", minWidth = 75),
      Focus = reactable::colDef(name = "Study focus", minWidth = 75, align = "right"),
      ARDS_definition = reactable::colDef(name = "ARDS definition", minWidth = 75, align = "right"),
      ARDS_pts = reactable::colDef(name = "Total ARDS patients", minWidth = 75, align = "right")
    )
  )

  return(table)
}

#' @title Create a table of study methods
#' @description Create an interactive html table of study methods
#' @param data -- Input in the format of `data_study`
#' @return An html table
#' @examples
#' \dontrun{
#' if(interactive()){
#'  methods_table(data_study)
#'  }
#' }
#' @export
#' @rdname methods_table
#'
#' @import dplyr
#' @import reactable
#' @importFrom reactablefmtr espn cell_style pill_buttons icon_assign

methods_table <- function(data) {

  ## Subset data

  data <- data |>
    dplyr::select(c(.data$First_author, .data$PMID, .data$Year, .data$Method, .data$Technology, .data$Tissue, .data$Cell))

  ## Themed html table

  table <- reactable::reactable(
    data,
    theme = reactablefmtr::espn(),
    fullWidth = TRUE,
    style = list(fontFamily = "Work Sans, sans-serif", fontSize = "20px"),
    showSortIcon = TRUE,
    showSortable = TRUE,
    searchable = FALSE,
    pagination = FALSE,
    defaultSorted = list(First_author = "asc", Year = "desc"),
    highlight = TRUE,
    borderless = TRUE,
    striped = TRUE,
    outlined = TRUE,
    resizable = TRUE,
    wrap = FALSE,
    columns = list(
      First_author = reactable::colDef(name = "Study", sticky = "left", minWidth = 75),
      PMID = reactable::colDef(name = "PMID", sticky = "left", minWidth = 75),
      Year = reactable::colDef(name = "Year", sticky = "left", minWidth = 75),
      Method = reactable::colDef(name = "Methodology", minWidth = 75, align = "right"),
      Technology = reactable::colDef(name = "Technology", minWidth = 75, align = "right"),
      Tissue = reactable::colDef(name = "Tissue type", minWidth = 100, align = "right"),
      Cell = reactable::colDef(name = "Cell type", minWidth = 100, na = "-", align = "right")
    )
  )

  return(table)
}
