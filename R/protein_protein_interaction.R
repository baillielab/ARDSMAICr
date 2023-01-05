#' @title Initialise STRING db
#' @description Initialises STRING db using the Bioconducter package `STRINGdb` and saves to a temp
#'     directory
#' @param vers -- STRING dbVersion -- Default: '11.5'
#' @param score_threshold -- Score threshold --  Default: 200
#' @param input_directory -- Input directory (can save for later use) -- Default: ''
#' @return An S4 object
#' @details Creates a full network
#' @examples
#' \dontrun{
#' if(interactive()){
#'  initialise_string_db()
#'  }
#' }
#' @export
#' @rdname initialise_string_db
#'
#' @import STRINGdb

initialise_string_db <- function(vers = "11.5", score_threshold = 200, input_directory = "") {

  ## Warning over memory issues

  message("This reads > 2Gb into memory...")

  continue <- readline(prompt = "Press 1 to continue or 2 to exit")
  continue_value <- as.integer(continue)

  if (continue == 2) {
    stop("STRING db initialisation halted...")
  } else {

    ## Intialises STRING db for H.sapiens to a temp dir

    STRINGdb::STRINGdb$new(
      version = vers,
      species = 9606,
      score_threshold = score_threshold,
      network_type = "full",
      input_directory = input_directory
    )
  }
}

#' @title Map genes to STRING IDs
#' @description Maps input genes to STRING IDs
#' @param string_db -- Object output from `initialise_string_db`
#' @param data_genes -- Input in the format of `data_genes`
#' @param n_genes P-- Number of genes to map from highest ranked maic gene to n
#' @return A tibble
#' @examples
#' \dontrun{
#' if(interactive()){
#'  initialise_string_db()
#'  mapped_genes <- map_string_db(data_genes, n_genes = 100)
#'  }
#' }
#' @export
#' @rdname map_string_db
#' @seealso see the STRINGdb package for further functions to work with the output
#'
#' @import dplyr

map_string_db <- function(string_db, data_genes, n_genes = NULL) {

  ## Set behaviour for number of genes to map

  if (is.null(n_genes)) {

    ## Default -- all maic genes

    example_mapped <- string_db$map(data_genes, "gene", removeUnmappedRows = TRUE)

    } else {

    ## User input

    data_genes_slice <- data_genes |>
        dplyr::slice(1:n_genes)
    example_mapped <- string_db$map(data_genes_slice, "gene", removeUnmappedRows = TRUE)
  }

  return(example_mapped)
}

