#' @title Create a chord diagram of contributions by study
#' @description Creates an interactive html plot of contributions by study
#' @param data -- Input in the format of `data_contributions`
#' @return An html chord plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  contributions_chord_bystudy(data_genes)
#'  }
#' }
#' @export
#' @rdname contributions_chord_bystudy
#'
#' @import dplyr
#' @import janitor
#' @import tibble
#' @import chorddiag

contributions_chord_bystudy <- function(data) {

  ## Subset and tidy data for plot

  circle <- data |>
    dplyr::select(c(.data$study, .data$contribution))

  by_study <- as.data.frame(t(circle)) |>
    janitor::row_to_names(row_number = 1) |>
    tibble::remove_rownames() |>
    dplyr::mutate_if(is.character, as.numeric)

  by_study <- as.matrix(by_study)

  row.names(by_study) <- ""

  ## Set colour scheme

  colours <- c(
    "white", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
           "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"
  )

  ## Html chord diagram

  plot <- chorddiag::chorddiag(
    by_study,
    type = "bipartite",
    groupColors = colours,
    showTicks = F,
    groupnameFontsize = 14,
    groupnamePadding = 15,
    margin = 200
  )

  return(plot)
}

#' @title Create a chord diagram of contributions by method type
#' @description Creates an interactive html plot of contributions by method type
#' @param data_contributions -- Input in the format of `data_contributions`
#' @param data_study -- Input in the format of `data_study`
#' @return An html chord plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  contributions_chord_bymethod(data_contributions, data_study)
#'  }
#' }
#' @export
#' @rdname contributions_chord_bymethod
#'
#' @import dplyr
#' @import forcats
#' @import janitor
#' @import tibble
#' @import chorddiag

contributions_chord_bymethod <- function(data_contributions, data_study) {

  ## Subset data

  data_study_grouped <- data_study |>
    dplyr::select(c(.data$uID, .data$Method))

  ## Rename column for join

  data_contributions_renamed <- data_contributions |>
    dplyr::rename(uID = "study")

  ## Join datasets by uID

  data_combined <- dplyr::left_join(data_contributions_renamed, data_study_grouped, by = "uID")

  ## Calculate contributions by method type

  by_type <- data_combined |>
    dplyr::mutate(Method = forcats::as_factor(.data$Method)) |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(total_contribution = sum(.data$contribution)) |>
    dplyr::mutate(total_contribution = .data$total_contribution * 100) |>
    dplyr::mutate_if(is.numeric, round, 3)

  ## Tidy data for plot

  by_type <- as.data.frame(t(by_type)) |>
    janitor::row_to_names(row_number = 1) |>
    tibble::remove_rownames() |>
    dplyr::mutate_if(is.character, as.numeric)

  by_type <- as.matrix(by_type)

  row.names(by_type) <- ""

  ## Html chord diagram

  plot <- chorddiag::chorddiag(
    by_type,
    type = "bipartite",
    groupColors = c("#bdbdbd", "#1b9e77", "#d95f02", "#7570b3"),
    showTicks = F,
    groupnameFontsize = 12,
    groupnamePadding = 10,
    margin = 120
  )

  return(plot)
}
