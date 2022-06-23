#' Generate pathfindR-like bubble plot from pathfindR results
#'
#' One odd choice in pathfindR is that the (very helpful) bubble plots it
#' generates get saved as an image file to your computer, but are otherwise
#' inaccessible to your R session. This function makes a similar bubble plot,
#' though with a bit more flexibility in how it's plotted and using a viridis
#' colour scale instead of the hard-to-read grey-red default.
#'
#' @param pathRes A dataframe, the output of dge_RunPathfinder() or, more
#'   likely, dge_CleanPathfinder().
#' @param nRes Numeric; number of pathways you would like to output.
#' @param xlim Length-2 numeric vector; limits of the x-axis (the fold-change;
#'   *technically* the y axis but we've flipped it here and so I'm naming as the
#'   graph outputs rather than behind the scenes...)
#' @param ylim Length-2 numeric vector; limits of the y-axis (not usually set;
#'   see nRes instead)
#' @param ... Additional args to pass to scale_color_viridis_c(), such as limits
#'   to set the bounds of the colour scale
#'
#' @return A ggplot object showing the fold enrichment, # genes, and p values of
#'   the top nRes pathways identified by pathfindeR.
#'
#' @examples
#'
#' pathRes <- system.file("extdata", "pathRes.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#'     # Raw output from dge_RunPathfinder:
#'     bubblePlot(pathRes)
#'
#' @importFrom rlang .data
#'
#' @export
#'
bubblePlot <- function(pathRes, nRes = 10, xlim = NULL, ylim = NULL, ...){
  if(!"Term_Description" %in% colnames(pathRes)){
    stop("Please ensure pathRes is a single dataframe, the output of
         dge_RunPathfinder. Vectorize as needed if your output is a list with
         multiple dataframes.")
  }
  pathRes |>
    dplyr::select(.data$Term_Description,
                  `Fold Enrichment` = "Fold_Enrichment",
                  .data$lowest_p,
                  .data$Up_regulated,
                  .data$Down_regulated) |>
    dplyr::slice_head(n = nRes) |>
    dplyr::mutate(Genes = paste(.data$Up_regulated, .data$Down_regulated, sep = ", "),
           `-Log10 P value` = -log10(.data$`lowest_p`)) |>
    dplyr::mutate(Genes = stringr::str_remove(.data$Genes, "^, ")) |>
    dplyr::mutate(`# Genes` = stringr::str_count(.data$Genes, ", ")) |>
    dplyr::mutate(Term_Description = forcats::fct_reorder(.data$Term_Description, .data$lowest_p, .desc = TRUE)) |>
    dplyr::select(-.data$Up_regulated, -.data$Down_regulated) |>
    ggplot2::ggplot(ggplot2::aes(.data$Term_Description, .data$`Fold Enrichment`, color = .data$`-Log10 P value`, size = .data$`# Genes`)) +
    ggplot2::geom_point() +
    ggplot2::coord_flip(xlim = ylim, ylim = xlim) +
    ggplot2::labs(x = ggplot2::element_blank()) +
    ggplot2::scale_color_viridis_c(...) +
    ggpubr::theme_pubr() +
    ggplot2::theme(legend.position = "right")
}
