#' @title Calculate the maic score inflection point
#' @description Calculates the maic score and gene number at which the inflection point in the curve
#'     of genes ranked by maic score occurs
#' @param data -- Input in the format of `data_genes`
#' @return A list of values
#' @examples
#' \dontrun{
#' if(interactive()){
#'  inflection_point(data_genes)
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

#' @title Create a plot of the inflection point
#' @description Create a plot of genes ranked by maic score with the inflection point annotated and
#'     point density of genes illustrated
#' @param data -- Input in the format of `data_genes`
#' @return A point density plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  inflection_point_plot(data_genes)
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

inflection_point_plot <- function(data) {

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

  ## Find the total number of genes

  total_gene_n <- nrow(data)

  ## Find the minimum maic score

  min_maic <- min(data$maic_score)

  ## Find the maximum maic score

  max_maic <- max(data$maic_score)

  ## Point density plot

  rank_plot <- ggplot2::ggplot(
    rank_data,
    ggplot2::aes(x = rank_data$rowname, y = rank_data$maic_score)
  ) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::geom_hline(yintercept = inflection_maic_score, linetype = "dashed", color = "steelblue") +
    ggplot2::geom_vline(xintercept = inflection_gene_n, linetype = "dashed", color = "steelblue") +
    viridis::scale_color_viridis("Point density") +
    ggplot2::scale_x_continuous(limits = c(0, total_gene_n), breaks = c(1, inflection_gene_n, seq(1000, total_gene_n, 500)), expand = c(0.05, 0)) +
    ggplot2::scale_y_continuous(limits = c(min_maic-0.2, max_maic+0.2), breaks = seq(min_maic, max_maic+0.2, 0.2), expand = c(0, 0)) +
    ggplot2::labs(x = "Gene rank", y = "MAIC score") +
    ggpubr::theme_pubr(legend = "right")

  return(rank_plot)
}

