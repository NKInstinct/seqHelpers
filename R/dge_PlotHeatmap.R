#' Plot CPM heatmap with sensible defaults
#'
#' Useful for getting a quick heatmap out of the output of dge_GetMatrix
#'
#' @param mat A matrix, likely the output of dge_GetMatrix() (not neccesarily,
#'   though).
#' @param scale Arg for ggheatmap, defaulting to "row" (gene-wise)
#' @param cluster_rows Arg for ggheatmap, defaulting to TRUE (yes, cluster genes)
#' @param show_cluster_rows Arg for ggheatmap, defaulting to FALSE (no, don't print the dendrogram)
#' @param ... Additional args to pass to ggheatmap
#'
#' @return A ggheatmap object
#' @export
#'

dge_PlotHeatmap <- function(mat, scale = "row", cluster_rows = TRUE, show_cluster_rows = FALSE, ...){
  gg <- ggheatmap::ggheatmap(mat, scale = scale, cluster_rows = cluster_rows, show_cluster_rows = show_cluster_rows, ...)
  return(gg)
}

#' @examples
#' dgeRes <- system.file("extdata", "dgeRes_full.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#' mat <- dge_GetMatrix(dgeRes$DGEList)
#'
#' dge_PlotHeatmap(mat)
#'
#'
#'
#'
#'
