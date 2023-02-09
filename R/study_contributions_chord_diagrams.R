#' @title Contributions chord by study
#' @description Creates an interactive html chord diagram of contributions to MAIC grouped by study.
#' @param data Data frame in the format of `ARDSMAICR::data_contributions`
#' @return An html chord diagram
#' @details
#' Input columns for `data_contributions` should be:
#' * `study` - Study ID - chr
#' * `raw_information` - Raw information. Sum of gene scores for each study - dbl
#' * `information` - Information relative to the sum of all gene scores - dbl
#' * `raw_contribution` - Raw contribution. Sum of gene scores for each study only where a score
#'     contributes to the overall MAIC score - dbl
#' * `contribution` - Contribution relative to the sum of all such scores - dbl
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  contributions_chord_bystudy(ARDSMAICR::data_contributions)
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
    dplyr::select(c(.data$study, .data$contribution)) |>
    dplyr::mutate(study = stringr::str_replace(.data$study, "(\\d{4})(\\d+)([^\\d]+)*", "\\1\\3")) |>
    dplyr::mutate(contribution = .data$contribution*100) |>
    dplyr::mutate(contribution = round(.data$contribution, 2))

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

#' @title Contributions chord by method
#' @description Creates an interactive html chord diagram of contributions to MAIC grouped by
#'     method.
#' @param data_contributions Data frame in the format of `ARDSMAICR::data_contributions`
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @return An html chord diagram
#' @details
#' Input columns for `data_contributions` should be:
#' * `study` - Study ID - chr
#' * `raw_information` - Raw information. Sum of gene scores for each study - dbl
#' * `information` - Information relative to the sum of all gene scores - dbl
#' * `raw_contribution` - Raw contribution. Sum of gene scores for each study only where a score
#'     contributes to the overall MAIC score - dbl
#' * `contribution` - Contribution relative to the sum of all such scores - dbl
#'
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
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  contributions_chord_bymethod(ARDSMAICR::data_contributions, ARDSMAICR::data_study)
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
