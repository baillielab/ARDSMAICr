#' @title Calculate study contributions
#' @description Calculates the raw and relative information an contribution of studies to MAIC.
#' @param data Data frame in the format of `ARDSMAICR::data_genes`
#' @return A tibble
#' @details
#' Input columns for `data_genes` should be (this is the standard output of the MAIC algorithm):
#' * `gene` - HGNC gene symbol - chr
#' * 1...
#' * `uID` - Study unique identifier. Column contains study specific gene score - dbl
#' * n...
#' * `maic_score` - MAIC score for gene - dbl
#' * `contributors` - Studies contributing to MAIC score by method - chr
#'
#' Definitions:
#' * `raw_information` = sum of gene scores for each study
#' * `relative_information` = `raw_information`/sum of all gene scores in MAIC
#' * `raw_contribution` = sum of gene scores for each study which contribute to the MAIC score
#' * `relative_contribution` = `raw_contribution`/sum of all gene scores in MAIC contributing to the
#' MAIC score
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  contributions_calculation(ARDSMAICR::data_genes)
#'  }
#' }
#' @export
#' @rdname contributions_calculation
#'
#' @import dplyr
#' @import tidyselect
#' @import tibble
#' @import tidyr
#' @import stringr

contributions_calculation <- function(data) {

  ## calculate raw information (sum of all gene scores for each study)

  raw_information <- data |>
    dplyr::select(c(tidyselect::where(is.double))) |>
    dplyr::select(-.data$maic_score) |>
    dplyr::summarise(across(tidyselect::where(is.numeric), sum)) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      tidyselect::everything(),
      names_to = "study",
      values_to = "raw_information"
    )

  ## calculate relative information (raw_information/sum of all gene scores in MAIC)

  relative_information <- raw_information |>
    dplyr::mutate(information = .data$raw_information / sum(.data$raw_information))

  ## calculate raw contribution (sum of all gene scores for each study which contributes to the MAIC score)

  raw_contribution <- data |>
    tidyr::separate(.data$contributors, c("study_1", "study_2", "study_3", "study_4", "study_5"), "\\,", fill = "right") |>
    dplyr::mutate_at(c("study_1", "study_2", "study_3", "study_4", "study_5"), ~ stringr::str_remove(.x, "^(.*?:)")) |>
    dplyr::mutate_at(c("study_1", "study_2", "study_3", "study_4", "study_5"), ~ stringr::str_trim(.x)) |>
    dplyr::mutate_at(c("study_1", "study_2", "study_3", "study_4", "study_5"), ~ tidyr::replace_na(., "")) |>
    dplyr::select(-.data$maic_score) |>
    dplyr::rowwise() |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ ifelse((
      dplyr::cur_column() %in% study_1 | dplyr::cur_column() %in% study_2 | dplyr::cur_column() %in% study_3 | dplyr::cur_column() %in% study_4 | dplyr::cur_column() %in% study_5),
      ., NA))) |>
    dplyr::ungroup() |>
    dplyr::summarise(across(tidyselect::where(is.numeric), sum, na.rm = TRUE)) |>
    tidyr::pivot_longer(
      tidyselect::everything(),
      names_to = "study",
      values_to = "raw_contribution"
    )

  ## calculate relative contribution (raw_contribution/sum of all gene scores contributing to the MAIC score)

  relative_contribution <- raw_contribution |>
    dplyr::mutate(contribution = .data$raw_contribution / sum(.data$raw_contribution))

  ## join information and contribution into single tibble

  contributions <- dplyr::left_join(relative_information, relative_contribution, by = "study") |>
    dplyr::arrange(dplyr::desc(.data$contribution))

  return(contributions)
}