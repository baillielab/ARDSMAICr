#' @title Count genes
#' @description Returns the total number of unique genes or SNPs included in MAIC.
#' @param data Data frame in the format of `ARDSMAICR::data_genes`
#' @return An integer
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
#'  count_genes(ARDSMAICR::data_genes)
#'  }
#' }
#' @export
#' @rdname count_genes

count_genes <- function(data) {

  ## converts to lowercase to maintain case sensitivity

  genes_unique <- unique(tolower(data$gene))

  length(genes_unique)
}





#' @title Count gene lists
#' @description Returns the total number of gene lists included in MAIC.
#' @param data Data frame in the format of `ARDSMAICR::data_study`
#' @return An integer
#' @details
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
#'  count_lists(ARDSMAICR::data_study)
#'  }
#' }
#' @export
#' @rdname count_lists
#'
#' @import dplyr

count_lists <- function(data) {

  ## Check there are no duplicate entries

  n_unique_id <- unique(data$uID)
  n_unique_id <- length(n_unique_id)
  n_list <- nrow(data)

  if (n_list > n_unique_id) {

    stop("Duplicate lists detected...")

  } else {

    ## Filter only available gene lists (included in MAIC)

    data <- data |>
      dplyr::mutate(List_available = as.logical(.data$List_available)) |>
      dplyr::filter(.data$List_available == TRUE)

    nrow(data)
  }
}





#' @title Count genes by gene list
#' @description Returns the number of unique genes or SNPs identified in each gene list (the most
#'     concise categorisation of lists).
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
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  print(genes_per_list(ARDSMAICR::data_genes), n = 50)
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





#' @title Count genes by method type
#' @description Returns the number of unique genes or SNPs identified by each method (the most
#'     concise categorisation of lists).
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @return A tibble
#' @details
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
#'  genes_per_method(ARDSMAICR::data_study, ARDSMAICR::data_genes)
#'  }
#' }
#' @export
#' @rdname genes_per_method
#'
#' @import dplyr
#' @import tidyselect
#' @import tidyr

genes_per_method <- function(data_study, data_genes) {

  ## Check for missing values in the methods column

  if (sum(is.na(data_study$Method)) > 0) {

    stop("Studies without a method detected...")

  } else {

    methods <- data_study |>
      dplyr::filter(.data$List_available == TRUE) |>
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
}





#' @title genes_per_category
#' @description Returns the number of unique genes or SNPs identified by each category defined in
#'     MAIC.
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @return A tibble
#' @details
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
#'  genes_per_category(ARDSMAICR::data_study, ARDSMAICR::data_genes)
#'  }
#' }
#' @rdname genes_per_category
#' @export
#'
#' @importFrom dplyr left_join filter select rename mutate across summarise arrange group_by ungroup
#' @importFrom tidyselect where everything
#' @importFrom tidyr pivot_longer

genes_per_category <- function(data_study, data_genes) {

  data_meta <- dplyr::left_join(data_study, categories, by = "uID")

  if (sum(is.na(data_meta$Method)) > 0) {

    stop("Studies without a method detected...")

  } else if (sum(is.na(data_meta$Category)) > 0) {

    stop("Studies without a category detected...")

  } else {

    categories_df <- data_meta |>
      dplyr::filter(.data$List_available == TRUE) |>
      dplyr::select(c(.data$uID, .data$Method, .data$Category)) |>
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

    n_gene_categories <- dplyr::left_join(data_gene_count, categories_df, by = "list")

    result <- n_gene_categories |>
      dplyr::group_by(.data$Category, .data$Method) |>
      dplyr::summarise(n_genes = sum(.data$n_genes)) |>
      dplyr::ungroup()

    return(result)
  }
}





#' @title Count gene lists by gene
#' @description For each unique gene or SNP returns the number of lists that gene or SNP is
#'     found in.
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
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ARDSMAICR::data_genes |>
#'    lists_per_gene() |>
#'    dplyr::filter(.data$gene == "IL6")
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





#' @title Count gene lists by method
#' @description Returns the total number of gene lists included in MAIC grouped by
#'     study method (the most concise categorisation of lists).
#' @param data Data frame in the format of `ARDSMAICR::data_study`
#' @return A tibble
#' @details
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
#'  lists_per_method(ARDSMAICR::data_study)
#'  }
#' }
#' @export
#' @rdname lists_per_method
#'
#' @import dplyr
#' @importFrom rlang .data

lists_per_method <- function(data) {

  ## Check for missing values in the methods column

  if (sum(is.na(data$Method)) > 0) {

    stop("Studies without a method detected...")

  } else {

    ## Filter available lists, group by method, count lists, and calculate %

    data |>
      dplyr::filter(.data$List_available == TRUE) |>
      dplyr::group_by(.data$Method) |>
      dplyr::summarise(count = dplyr::n()) |>
      dplyr::mutate(percentage = .data$count / sum(.data$count))
  }
}





#' @title Count gene lists by method by gene
#' @description Returns the number of lists each unique gene or SNP is found in stratified by
#'     method (the most concise categorisation of lists).
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @return A tibble
#' @details
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
#'    lists_per_method_per_gene(ARDSMAICR::data_study, ARDSMAICR::data_genes) |>
#'      dplyr::filter(.data$gene == "IL6")
#'  }
#' }
#' @export
#' @rdname lists_per_method_per_gene
#'
#' @import dplyr
#' @import tidyselect
#' @import tidyr

lists_per_method_per_gene <- function(data_study, data_genes) {

  ## Check for missing values in the methods column

  if (sum(is.na(data_study$Method)) > 0) {

    stop("Studies without a method detected...")

  } else {

  methods <- data_study |>
    dplyr::filter(.data$List_available == TRUE) |>
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
    dplyr::summarise(n_lists = sum(.data$identified)) |>
    dplyr::ungroup()
  }
}





#' @title Count gene lists by category
#' @description Returns the total number of gene lists included in MAIC grouped by
#'     study category.
#' @param data Data frame in the format of `ARDSMAICR::data_study`
#' @return A tibble
#' @details
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
#'  lists_per_category(ARDSMAICR::data_study)
#'  }
#' }
#' @rdname lists_per_category
#' @export
#'
#' @importFrom dplyr left_join filter group_by summarise n ungroup mutate

lists_per_category <- function(data) {

  data_meta <- dplyr::left_join(data, categories, by = "uID")

  if (sum(is.na(data_meta$Method)) > 0) {

    stop("Studies without a method detected...")

  } else if (sum(is.na(data_meta$Category)) > 0) {

    stop("Studies without a category detected...")

  } else {

    result <- data_meta |>
      dplyr::filter(.data$List_available == TRUE) |>
      dplyr::group_by(.data$Category, .data$Method) |>
      dplyr::summarise(count = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(percentage = 100*(.data$count / sum(.data$count)))

    return(result)

  }
}





#' @title Count gene lists by category by gene
#' @description Returns the number of lists each unique gene or SNP is found in stratified by
#'    category.
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @return A tibble
#' @details
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
#'  lists_per_category_per_gene(ARDSMAICR::data_study, ARDSMAICR::data_genes) |>
#'    dplyr::filter(gene == "IL6")
#'  }
#' }
#' @rdname lists_per_category_per_gene
#' @export
#'
#' @importFrom dplyr left_join filter select mutate across group_by summarise ungroup
#' @importFrom tidyselect where
#' @importFrom tidyr pivot_longer

lists_per_category_per_gene <- function(data_study, data_genes) {

  data_meta <- dplyr::left_join(data_study, categories, by = "uID")

  ## Check for missing values in the methods column

  if (sum(is.na(data_study$Method)) > 0) {

    stop("Studies without a method detected...")

  } else if (sum(is.na(data_meta$Category)) > 0) {

    stop("Studies without a category detected...")

  } else {

    categories_df <- data_meta |>
      dplyr::filter(.data$List_available == TRUE) |>
      dplyr::select(c(.data$uID, .data$Method, .data$Category))

    data_categories_count <- data_genes |>
      dplyr::select(-c(.data$maic_score, .data$contributors)) |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
      tidyr::pivot_longer(cols = tidyselect::where(is.logical),
                          names_to = "uID",
                          values_to = "identified"
      )

    n_categories_gene <- dplyr::left_join(data_categories_count, categories_df, by = "uID", copy = TRUE)

    result <- n_categories_gene |>
      dplyr::group_by(.data$gene, .data$Category, .data$Method) |>
      dplyr::summarise(n_lists = sum(.data$identified)) |>
      dplyr::ungroup()

    return(result)
  }
}





#' @title Count methods by gene
#' @description Returns the number of methods (the most concise categorisation of lists) supporting
#'     each unique gene or SNP.
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
#' @md
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ARDSMAICR::data_genes |>
#'    methods_per_gene() |>
#'    dplyr::filter(.data$gene == "IL6")
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

  warning("Only GWAS, TRANSCRIPTOMICS, & PROTEOMICS supported as methods...")

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





#' @title Count categories per gene
#' @description Returns the number of categories supporting each unique gene or SNP.
#' @param data_study Data frame in the format of `ARDSMAICR::data_study`
#' @param data_genes Data frame in the format of `ARDSMAICR::data_genes`
#' @return A tibble
#' @details
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
#'  categories_per_gene(ARDSMAICR::data_study, ARDSMAICR::data_genes) |>
#'    dplyr::filter(gene == "IL6")
#'  }
#' }
#' @rdname categories_per_gene
#' @export
#'
#' @importFrom dplyr left_join filter select mutate across group_by summarise ungroup summarize
#' @importFrom tidyselect where
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_glue str_count
categories_per_gene <- function(data_study, data_genes) {

  data_meta <- dplyr::left_join(data_study, categories, by = "uID")

  ## Check for missing values in the methods column

  if (sum(is.na(data_study$Method)) > 0) {

    stop("Studies without a method detected...")

  } else if (sum(is.na(data_meta$Category)) > 0) {

    stop("Studies without a category detected...")

  } else {

    categories_df <- data_meta |>
      dplyr::filter(.data$List_available == TRUE) |>
      dplyr::select(c(.data$uID, .data$Method, .data$Category))

    data_categories_count <- data_genes |>
      dplyr::select(-c(.data$maic_score, .data$contributors)) |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.double), ~ as.logical(., na.rm = TRUE))) |>
      tidyr::pivot_longer(cols = tidyselect::where(is.logical),
                          names_to = "uID",
                          values_to = "identified"
      )

    n_categories_gene <- dplyr::left_join(data_categories_count, categories_df, by = "uID", copy = TRUE)

    result_long <- n_categories_gene |>
      dplyr::group_by(.data$gene, .data$Category, .data$Method) |>
      dplyr::summarise(n_lists = sum(.data$identified)) |>
      dplyr::ungroup()

    result <- result_long |>
      dplyr::filter(.data$n_lists >= 1) |>
      dplyr::mutate(Categories = stringr::str_glue("{Method}:{Category}")) |>
      dplyr::select(-c(.data$Category, .data$Method)) |>
      dplyr::group_by(.data$gene) |>
      dplyr::summarize(support = paste(.data$Categories, collapse = ",\ ")) |>
      dplyr::ungroup() |>
      dplyr::mutate(n_categories = 1 + stringr::str_count(.data$support, ","))

    return(result)
  }
}






