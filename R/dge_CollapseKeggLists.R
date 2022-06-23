#' Collapse a set of KeggGeneLists into one list of vectors
#'
#' The idea here is that if you run a DGE experiment with multiple comparisons,
#' you'll probably be vectorizing your whole analysis. However, at the end, when
#' you want to show heatmaps of significant genes in pathways, it'll be helpful
#' to take your list of significant KEGG terms (one set per comparison) and make
#' a single, main list that contains all significant genes for all signficiant
#' pathways, so that your heatmap (which will contain all samples) is
#' representative. This function does exactly that.
#'
#' @param keggResList List of results from running multiple
#'   dge_GetKeggGeneLists() commands, likely from vectorizing it over multiple
#'   dge_RunPathview requests.
#'
#' @return A single list of named vectors, like the output of
#'   dge_GetKeggGeneLists: one vector per pathway, contains all sig genes in the
#'   pathway.
#' @export
#'
dge_CollapseKeggLists <- function(keggResList){
  if(methods::is(keggResList[[1]], "data.frame")){
    stop("Unexpected dataframe input. Please ensure you are running this
         function on a list of outputs from dge_GetKeggGeneLists(), not on a
         pathfinder results table directly.")
  }
  res <- purrr::reduce(keggResList, \(a, b){
    keys <- unique(c(names(a), names(b)))
    purrr::map2(a[keys], b[keys], c) |>
      purrr::set_names(keys)
  })

  res <- res |>
    purrr::map(unique) |>
    purrr::map(stringr::str_subset, pattern = ".+")

  return(res)
}

#' @examples
#' pathRes_full<- system.file("extdata", "pathRes_full.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#' # In this case, pathRes_full contains the results of running
#' # dge_RunPathfinder() on all comparisons in pathRes_full.Rds (which we can't
#' # do here because that'll write reports into the package).
#'
#' multiCompLists <- purrr::map(pathRes_full, dge_GetKeggGeneLists)
#'
#' # To combine all the genes in each list into one, we just run:
#'
#' combinedList <- dge_CollapseKeggLists(multiCompLists)
