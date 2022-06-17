#' Get a list of all significantly altered genes in all significant pathways
#'
#' Returns one list per pathway, from the output of dge_RunPathfinder().
#'
#' @param pathRes Dataframe, the output of dge_RunPathfinder().
#'
#' @return A named list of character vectors containing gene names, with each
#'   vector named by its KEGG pathway.
#' @importFrom rlang .data
#' @export
#'
dge_GetKeggGeneLists <- function(pathRes){
  df <- dplyr::select(pathRes, Term = "Term_Description", Up = "Up_regulated", Down = "Down_regulated") |>
    dplyr::mutate(Genes = paste(.data$Up, .data$Down, sep = ", "),
           Genes = stringr::str_remove(.data$Genes, "^, ")) |>
    dplyr::select(.data$Term, .data$Genes)

  geneList <- lapply(with(df, split(Genes, Term)), as.list) |>
    unlist() |>
    purrr::map(stringr::str_split, pattern = ", ") |>
    purrr::map(unlist)

  return(geneList)
}

#' @examples
#'
#' pathRes <- system.file("extdata", "pathRes.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#' dge_GetKeggGeneLists(pathRes)
