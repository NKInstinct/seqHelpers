#' Clean up the output of dge_RunPathfinder to make nicer bubble plots
#'
#' Small cleaning function to do things like wrap long pathway names, and filter
#' on fold enrichment scores and gene counts to get rid of irrelevant pathways.
#'
#' @param pathRes Dataframe, output of dge_RunPathfinder().
#' @param wrap Numeric, max number of characters in a pathway name before
#'   wrapping to a new line.
#' @param foldChange Numeric, keep only pathways with a fold enrichment **over**
#'   this number.
#' @param geneCount Numeric, keep only pathways with a gene count **over** this
#'   number.
#'
#' @return Dataframe with the above cleaning steps applied.
#'
#' @importFrom rlang .data
#' @export
#'

dge_CleanPathfinder <- function(pathRes, wrap = 35, foldChange = 1, geneCount = 20){
  if(!"Term_Description" %in% colnames(pathRes)){
    stop("Please ensure pathRes is a single dataframe, the output of
         dge_RunPathfinder. Vectorize as needed if your output is a list with
         multiple dataframes.")
  }
  df <- pathRes |>
    dplyr::mutate(Term_Description = stringr::str_wrap(.data$Term_Description, width = wrap)) |>
    dplyr::mutate(GeneCount = paste(.data$Up_regulated, .data$Down_regulated, sep = ", "),
           GeneCount = stringr::str_remove(.data$GeneCount, "^, "),
           GeneCount = stringr::str_count(.data$GeneCount, pattern = stringr::boundary("word"))) |>
    dplyr::filter(.data$Fold_Enrichment > foldChange,
                  .data$GeneCount > geneCount) |>
    dplyr::select(-.data$GeneCount)
}

#' @examples
#'
#' pathRes <- system.file("extdata", "pathRes.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#' pathClean <- dge_CleanPathfinder(pathRes, geneCount = 2)
