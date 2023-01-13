#' @title MAIC overlap table
#' @description Creates an interactive html table of genes found in common between MAIC and a
#'     and a user supplied list of genes.
#' @param data_genes Data frame in the format of `data_genes`
#' @param data_alternative Data frame where the only column is gene symbols and is titled gene,
#'     unless is `data_biolitmine`
#' @param biolitmine Boolean -- Default = TRUE -- if comparing BioLitMine result using
#'     `data_biolitmine`
#' @return An html table
#' @details
#' Primary use is to compare ARDS MAIC with the results of a BioLitMine search for the MeSH
#' term "Respiratory Distress Syndrome, Adult". A secondary use case is to compare ARDS MAIC
#' with MAIC analyses performed in different conditions e.g., our previous MAIC of SARS-CoV-2
#' whole-genome studies.
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
#'  overlap_table(data_genes, data_covidmaic, biolitmine = FALSE)
#'  }
#' }
#' @export
#' @rdname overlap_table
#'
#' @import dplyr
#' @import tibble
#' @import reactable
#' @importFrom reactablefmtr espn

overlap_table <- function(data_genes, data_alternative, biolitmine = TRUE) {

  ## Extract maic gene names

  maic_gene_names <- data_genes$gene

  ## Set biolitmine behaviour

  if (biolitmine == TRUE) {

    ## Check biolitmine format

    col_names <- colnames(data_alternative)

    check_names <- c("Gene", 'Last_10_PMIDS', "Publication_count")

    diff <- dplyr::setdiff(check_names, col_names)

    if (length(diff) >= 1) {

      stop("Data not in BiolLitMine format...")

    }

    ## Find common genes between maic and biolitmine

    common_genes <- maic_gene_names[(maic_gene_names %in% data_alternative$Gene)]

    ## Filter data_biolitmine for common genes

    common_alternative <- data_alternative |>
      dplyr::filter(.data$Gene %in% common_genes) |>
      dplyr::arrange(desc(.data$Publication_count)) |>
      dplyr::rename(gene = .data$Gene)

    ## Filter data_genes for common genes and select gene name and maic ranking

    common_data_genes <- data_genes |>
      tibble::rowid_to_column(var = "maic_ranking") |>
      dplyr::filter(.data$gene %in% common_genes) |>
      dplyr::select(c(.data$gene, .data$maic_ranking))

    ## Join data frames

    common <- left_join(common_alternative, common_data_genes, by = "gene")

    ## Themed html table

    table <- reactable::reactable(
      common,
      theme = reactablefmtr::espn(),
      fullWidth = TRUE,
      style = list(fontFamily = "Work Sans, sans-serif", fontSize = "20px"),
      showSortIcon = TRUE,
      showSortable = TRUE,
      searchable = FALSE,
      pagination = FALSE,
      defaultSorted = list(Publication_count = "desc", gene = "asc"),
      highlight = TRUE,
      borderless = TRUE,
      striped = TRUE,
      outlined = TRUE,
      resizable = TRUE,
      wrap = FALSE,
      columns = list(
        gene = reactable::colDef(name = "Gene", sticky = "left", minWidth = 75),
        Last_10_PMIDS = reactable::colDef(name = "Last 10 publications", sticky = "left", minWidth = 125),
        Publication_count = reactable::colDef(name = "Publication count", sticky = "left", minWidth = 75),
        maic_ranking = reactable::colDef(name = "MAIC Rank", sticky = "left", minWidth = 75)
      )
    )
   } else {

    ## Find common genes between maic and alternative gene list

    common_genes <- maic_gene_names[(maic_gene_names %in% dplyr::pull(data_alternative[, 1]))]

    ## Filter data_genes for common genes and select gene name and maic ranking

    common <- data_genes |>
      tibble::rowid_to_column(var = "maic_ranking") |>
      dplyr::filter(.data$gene %in% common_genes) |>
      dplyr::select(c(.data$gene, .data$maic_ranking))

    table <- reactable::reactable(
      common,
      theme = reactablefmtr::espn(),
      fullWidth = TRUE,
      style = list(fontFamily = "Work Sans, sans-serif", fontSize = "20px"),
      showSortIcon = TRUE,
      showSortable = TRUE,
      searchable = FALSE,
      pagination = FALSE,
      defaultSorted = list(maic_ranking = "asc", gene = "asc"),
      highlight = TRUE,
      borderless = TRUE,
      striped = TRUE,
      outlined = TRUE,
      resizable = TRUE,
      wrap = FALSE,
      columns = list(
        gene = reactable::colDef(name = "Gene", sticky = "left", minWidth = 75),
        maic_ranking = reactable::colDef(name = "MAIC Rank", sticky = "left", minWidth = 75)
      )
    )
  }

  return(table)
}

#' @title MAIC unidentified table
#' @description Create an interactive html table of genes found in an alternative list but not in
#'     MAIC.
#' @param data_genes Data frame in the format of `data_genes`
#' @param data_alternative Data frame where the only column is gene symbols and titled gene,
#'     unless is `data_biolitmine`
#' @param biolitmine Boolean -- Default = TRUE -- if comparing BioLitMine result using
#'     `data_biolitmine`
#' @return An html table (biolitmine) or a tibble (other source)
#' @details
#' Primary use is to compare ARDS MAIC with the results of a BioLitMine search for the MeSH
#' term "Respiratory Distress Syndrome, Adult". A secondary use case is to compare ARDS MAIC
#' with MAIC analyses performed in different conditions e.g., our pervious MAIC of SARS-CoV-2
#' whole-genome studies.
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
#'  overlap_unidentified_table(data_genes, data_biolitmine, biolitmine = TRUE)
#'  }
#' }
#' @export
#' @rdname overlap_unidentified_table
#'
#' @import dplyr
#' @import reactable
#' @importFrom reactablefmtr espn

overlap_unidentified_table <- function(data_genes, data_alternative, biolitmine = TRUE) {

  ## Extract maic gene names

  maic_gene_names <- data_genes$gene

  ## Set biolitmine behaviour

  if (biolitmine == TRUE) {

    ## Check biolitmine format

    col_names <- colnames(data_alternative)

    check_names <- c("Gene", 'Last_10_PMIDS', "Publication_count")

    diff <- dplyr::setdiff(check_names, col_names)

    if (length(diff) >= 1) {

      stop("Data not in BiolLitMine format...")

    }

    ## Find genes present in biolitmine not identified by maic

    unidentified_genes <- data_alternative$Gene[!(data_alternative$Gene %in% maic_gene_names)]

    ## Filter data_alternative for unidentified genes

    unidentified_alternative <- data_alternative |>
      dplyr::filter(.data$Gene %in% unidentified_genes) %>%
      dplyr::arrange(desc(.data$Publication_count))

    ## Themed html table

    table <- reactable::reactable(
      unidentified_alternative,
      theme = reactablefmtr::espn(),
      fullWidth = TRUE,
      style = list(fontFamily = "Work Sans, sans-serif", fontSize = "20px"),
      showSortIcon = TRUE,
      showSortable = TRUE,
      searchable = FALSE,
      pagination = FALSE,
      defaultSorted = list(Publication_count = "desc", Gene = "asc"),
      highlight = TRUE,
      borderless = TRUE,
      striped = TRUE,
      outlined = TRUE,
      resizable = TRUE,
      wrap = FALSE,
      columns = list(
        Gene = reactable::colDef(name = "Gene", sticky = "left", minWidth = 75),
        Last_10_PMIDS = reactable::colDef(name = "Last 10 publications", sticky = "left", minWidth = 125),
        Publication_count = reactable::colDef(name = "Publication count", sticky = "left", minWidth = 75)
      )
    )

    return(table)

  } else {

    ## Rename data_alternative column 1 as gene

    data_alternative <- data_alternative |>
      dplyr::rename(Gene = 1)

    ## Find genes present in biolitmine not identified by maic

    unidentified_genes <- data_alternative$Gene[!(data_alternative$Gene %in% maic_gene_names)]

    ## As tibble

    unidentifed_alternative <- tibble::as_tibble(unidentified_genes) %>%
      dplyr::rename(gene = .data$value)

    return(unidentifed_alternative)
  }
}

#' @title MAIC overlap Venn diagram
#' @description Creates a Venn diagram of the overlap between MAIC genes and genes found in an
#'     alternative list.
#' @param data_genes Data frame in the format of `data_genes`
#' @param data_alternative Data frame where the only column is gene symbols and titled gene,
#'     unless is `data_biolitmine`
#' @param biolitmine Boolean -- Default = TRUE -- if comparing BioLitMine result using
#'     `data_biolitmine`
#' @return Venn diagram
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
#'  overlap_venn(data_genes, data_biolitmine, biolitmine = TRUE)
#'  }
#' }
#' @export
#' @rdname overlap_venn
#'
#' @import dplyr
#' @import ggvenn

overlap_venn <- function(data_genes, data_alternative, biolitmine = TRUE) {

  ## Extract maic gene names

  maic_genes <- data_genes |>
    dplyr::pull(.data$gene)

  ## Set biolitmine behaviour

  if (biolitmine == TRUE) {

    ## Check biolitmine format

    col_names <- colnames(data_alternative)

    check_names <- c("Gene", 'Last_10_PMIDS', "Publication_count")

    diff <- dplyr::setdiff(check_names, col_names)

    if (length(diff) >= 1) {

      stop("Data not in BiolLitMine format...")

    }

    ## Extract biolitmine gene names

    alternative_genes <- data_alternative |>
      dplyr::pull(.data$Gene)

    ## Create list object

    venn_list <- list(`MAIC` = maic_genes, `ARDS literature` = alternative_genes)

    ## Plot venn diagram

    venn <- ggvenn::ggvenn(
      venn_list,
      c("MAIC", "ARDS literature"),
      fill_color = c("#66c2a5", "#fc8d62"),
      fill_alpha = 0.3,
      stroke_alpha = 0.8,
      stroke_size = 0.5,
      text_size = 5,
      show_percentage = FALSE
    )

  } else {

    ## Extract alternative list gene names

    alternative_genes <- dplyr::pull(data_alternative[, 1])

    ## Create list object

    venn_list <- list(`MAIC` = maic_genes, `Alt list` = alternative_genes)

    ## Plot venn diagram

    venn <- ggvenn::ggvenn(
      venn_list,
      c("MAIC", "Alt list"),
      fill_color = c("#66c2a5", "#fc8d62"),
      fill_alpha = 0.3,
      stroke_alpha = 0.8,
      stroke_size = 0.5,
      text_size = 5,
      show_percentage = FALSE
    )
  }

  return(venn)
}

#' @title Count MAIC and alternative gene list overlap
#' @description Returns the total number of genes identified in the alternative list, the number of
#' these genes found in MAIC, and the percentage of alternative list genes found in MAIC.
#' @param data_genes Data frame in the format of `data_genes`
#' @param data_alternative Data frame where the first column is a gene list
#' @return A list
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
#'  overlap_count(data_genes, data_biolitmine)
#'  }
#' }
#' @export
#' @rdname overlap_count

overlap_count <- function(data_genes, data_alternative) {

  ## Count number of genes found by biolitmine

  alternative_len <- nrow(data_alternative)

  ## Extract maic gene names

  maic_gene_names <- data_genes$gene

  ## Find common genes between maic and biolitmine

  common_genes <- maic_gene_names[(maic_gene_names %in% dplyr::pull(data_alternative[, 1]))]

  ## Count

  common_genes_len <- length(common_genes)

  ## Percentage

  percent <- round((common_genes_len / alternative_len) * 100, digits = 2)

  ## Build list

  values <- list(alternative_len, common_genes_len, percent)

  names(values) <- c("alternative_n", "maic_overlap_n", "maic_overlap_percentage")

  return(values)
}
