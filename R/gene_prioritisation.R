#' @title Calculate the inflection point
#' @description Calculates the MAIC score and gene number at which the inflection point in the curve
#'     of genes ranked by MAIC score occurs.
#' @param data Data frame in the format of `ARDSMAICR::data_genes`
#' @return A list of values
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
#'  inflection_point(ARDSMAICR::data_genes)
#'  }
#' }
#' @export
#' @rdname inflection_point
#' @references \url{https://CRAN.R-project.org/package=inflection}
#'
#' @import tibble
#' @import dplyr
#' @import inflection

inflection_point <- function(data) {

  ## Tidy data and sort by maic score

  rank_data <- data |>
    tibble::rownames_to_column(var = "rowname") |>
    dplyr::mutate(rowname = as.numeric(.data$rowname)) |>
    dplyr::arrange(dplyr::desc(.data$maic_score))

  ## Calculate the inflection point using the unit invariant knee method implemented in the
  ## `inflection` package

  inflection_gene_n <- inflection::uik(rank_data$rowname, rank_data$maic_score)

  ## Find the maic score and gene number at the inflection point

  maic_score <- rank_data |>
    dplyr::filter(.data$rowname == inflection_gene_n) |>
    dplyr::select(.data$maic_score)

  inflection_maic_score <- maic_score$maic_score

  values <- list(inflection_maic_score, inflection_gene_n)

  names(values) <- c("maic_score", "gene_number")

  return(values)
}

#' @title Plot the inflection point
#' @description Creates a plot of genes ranked by MAIC score with the inflection point annotated and
#'     the point density of genes illustrated.
#' @param data Data frame in the format of `ARDSMAICR::data_genes`
#' @param first_break Integer value for first break on x-axis after inflection point -- Default = 1000
#' @param increment Integer value for size of increment for subsequent breaks on x-axis -- Default = 500
#' @return A point density plot
#' @details
#' Input columns for `data_genes` should be (this is the standard output of the MAIC algorithm):
#' * `gene` - HGNC gene symbol - chr
#' * 1...
#' * `uID` - Study unique identifier. Column contains study specific gene score - dbl
#' * n...
#' * `maic_score` - MAIC score for gene - dbl
#' * `contributors` - Studies contributing to MAIC score by method - chr
#' @examples
#' \dontrun{
#' if(interactive()){
#'  inflection_point_plot(ARDSMAICR::data_genes, first_break = 2000, incremnt = 1000)
#'  }
#' }
#' @export
#' @rdname inflection_point_plot
#'
#' @import tibble
#' @import dplyr
#' @import inflection
#' @import ggplot2
#' @import ggpointdensity
#' @import viridis
#' @import ggpubr

inflection_point_plot <- function(data, first_break = 1000, increment = 500) {

  ## Tidy data and sort by maic score

  rank_data <- data |>
    tibble::rownames_to_column(var = "rowname") |>
    dplyr::mutate(rowname = as.numeric(.data$rowname)) |>
    dplyr::arrange(dplyr::desc(.data$maic_score)) |>
    dplyr::mutate(maic_score = round(.data$maic_score, digits = 3))


  ## Calculate the inflection point using the unit invariant knee method implemented in the
  ## `inflection` package

  inflection_gene_n <- inflection::uik(rank_data$rowname, rank_data$maic_score)

  ## Find the maic score and gene number at the inflection point

  maic_score <- rank_data |>
    dplyr::filter(.data$rowname == inflection_gene_n) |>
    dplyr::select(.data$maic_score) |>
    dplyr::mutate(maic_score = round(.data$maic_score, digits = 3))

  inflection_maic_score <- maic_score$maic_score

  ## Find the total number of genes

  total_gene_n <- nrow(data)

  ## Find the minimum maic score

  min_maic <- min(round(data$maic_score, digits = 3))

  ## Find the maximum maic score

  max_maic <- max(round(data$maic_score, digits = 3))

  ## Point density plot

  rank_plot <- ggplot2::ggplot(
    rank_data,
    ggplot2::aes(x = rank_data$rowname, y = rank_data$maic_score)
  ) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::geom_hline(yintercept = inflection_maic_score, linetype = "dashed", color = "steelblue") +
    ggplot2::geom_vline(xintercept = inflection_gene_n, linetype = "dashed", color = "steelblue") +
    viridis::scale_color_viridis("Point density") +
    ggplot2::scale_x_continuous(limits = c(0, total_gene_n), breaks = c(1, inflection_gene_n, seq(first_break, total_gene_n, increment)), expand = c(0.05, 0)) +
    ggplot2::scale_y_continuous(limits = c(min_maic-0.2, max_maic+0.2), breaks = seq(min_maic, max_maic+0.2, 0.2), expand = c(0, 0)) +
    ggplot2::labs(x = "Gene rank", y = "MAIC score") +
    ggpubr::theme_pubr(legend = "right")

  return(rank_plot)
}

