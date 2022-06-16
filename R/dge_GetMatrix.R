#' Get CPM-corrected count matrix out of a DGEList, then filter on genes of interest
#'
#' Get the DGEList out of a dge_OneFactor result, and return a matrix of
#' (selected) genes.
#'
#' @param DGEList A DGEList object, usually returned when you run
#'   dge_OneFactor() with retainDGEList = TRUE.
#' @param genes Character vector of gene names, in the same format as the
#'   rownames from your DGEList. If you're unsure what these will look like in
#'   your own data, run this function keeping genes = NULL and then peek
#'   (head()) at the output.
#' @param scale One of "none", "cols", or "rows"; should centering and scaling
#'   be applied to the data, and if so, should it be row-wise or column-wise?
#' @param colnames Optional vector of new colnames. This will be directly
#'   replaced 1:1, so it needs to be the same length as the number of columns in
#'   the output, and in the same order. Useful for when your sample names are
#'   gross and you want to replace them with something nicer for the resulting
#'   heatmap.
#'
#' @return A matrix with (selected) genes as rows, samples as columns, and
#'   values as CPM-adjusted counts with optional scaling.
#'
#' @examples
#'
#'   dgeRes <- system.file("extdata", "dgeRes_full.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#'   dge_GetMatrix(dgeRes$DGEList)
#'
#' @export
#'
dge_GetMatrix <- function(DGEList, genes = NULL, scale = "none", colnames = NULL){
  mat <- edgeR::cpm(DGEList)

  if(!is.null(genes)){
    mat <- mat[rownames(mat) %in% genes,]
  }
  if(scale == "cols"){
    mat <- scale(mat)
  }else if(scale == "rows"){
    mat <- mat |>
      t() |>
      scale() |>
      t()
  }

  if(!is.null(colnames)){
    colnames(mat) <- colnames
  }
  return(mat)
}
