#' Run pathfindR analysis on output from dge_OneFactor
#'
#' Essentially just a wrapper around run_pathfindR(), which unfortunately also
#' means it'll produce an HTML report and there's no way to stop it.
#'
#' @param dgeRes Dataframe; the output of running dge_OneFactor()
#' @param geneNames Column ID to choose the name of the column containing gene
#'   names (tidy evaluation).
#' @param species Character. Currently only "Mouse" is supported, but others can
#'   be added as needed.
#' @param reportFolder Path to folder where the pathfindR report will be saved.
#'
#' @return A dataframe showing all significantly altered KEGG pathways and the
#'   constituent genes that are up- and down-regulated.
#'
#' @examples
#'
#' dgeRes <- system.file("extdata", "dgeRes.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'  \dontrun{
#'     dge_RunPathfinder(dgeRes)
#'  }
#'
#' @importFrom rlang .data
#' @export
dge_RunPathfinder <- function(dgeRes, geneNames = .data$Newid, species = "Mouse", reportFolder = "pathfindR_Results"){
  if(species != "Mouse"){
    stop("Currently, only Mouse KEGG pathways are supported. Please submit an issue if you need other species")
  }
  geneNames <- rlang::enquo(geneNames)
  res <- dplyr::select(dgeRes, !!geneNames, .data$logFC, .data$PValue) |>
    pathfindR::run_pathfindR(convert2alias = FALSE,
                       gene_sets = "mmu_KEGG",
                       pin_name_path = "mmu_STRING",
                       list_active_snw_genes = FALSE,
                       output_dir = reportFolder)
  return(res)
}

