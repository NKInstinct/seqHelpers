#' Plot all significant genes in a heatmap, with GO-annotated clusters
#'
#' This function takes the output from dge_OneFactor, finds all significant
#' genes in all comparisons, and then runs GO term enrichment and annotates a
#' heatmap with cluster-specific GO enrichments.
#'
#' This strategy is designed to reduce the implicit bias GO analysis has if you
#' use "the whole genome" as your gene universe. Since the cells you are
#' comparing come from the same tissue, any false positive genes (i.e. 5-10% of
#' your significant list) are likely to be cell type/tissue relevant, and so are
#' likely to enrich for "relevant" looking GO terms. This approach instead uses
#' GO terms only to annotate clusters of similarly-expressed genes, using only
#' all significant genes as the gene universe.
#'
#' @param dgeRes List, output of dge_OneFactor. Note that this must contain a
#'   DGEList object!
#' @param genePVal PValue cutoff for identifying a gene as significant. Advise
#'   0.1 since the GO term enrichment is robust to false positives in this
#'   approach, and works better with a larger pool of genes to draw on.
#' @param geneID The column in dgeRes that contains the gene IDs. Usually Newid
#'   if name conversion was done and GeneID if not.
#' @param nClusters Number of clusters to use for unsupervised k-medoid
#'   clustering. Needs to be data-driven, so try something in the 5-10 range and
#'   then re-run if necessary to get the number that best fits your data. Can
#'   also map this if you like but the k-medoids algorithm is slow so it'll take
#'   a while.
#' @param convertMatrix Boolean, should the matrix gene IDs be converted? See
#'   dge_GetMatrix for more.
#' @param ... Additional arguments to pass to convertBiomart (i.e. convert_from
#'   and convert_to). See convertBiomart for more.
#'
#' @return A ggplot heatmap of all significant genes, clustered, and annotated
#'   with cluster-specific GO enrichments.
#'
#' @examples
#' dgeRes_full<- system.file("extdata", "dgeRes_full.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#' \dontrun{
#' dge_PlotGOHeatmap(dgeRes_full, genePVal = 0.1, geneID = Newid, nClusters = 2, convertMatrix = TRUE)
#' }
#'
#' @importFrom rlang .data
#' @import topGO
#' @export
dge_PlotGOHeatmap <- function(dgeRes, genePVal, geneID, nClusters, convertMatrix = FALSE, ...){
  if(is.null(dgeRes$DGEList)){
    stop("dgeRes must contain a DGEList entry. Please run dge_OneFactor with retainDGEList == TRUE")
  }

  geneID <- rlang::enquo(geneID)

  resList <- dgeRes[!names(dgeRes) == "DGEList"]

  sigMat <- collapseSigGenes(resList, genePVal, !!geneID) |>
    dge_GetMatrix(dgeRes$DGEList, scale = "rows", genes = _, convertMatrix = convertMatrix, ...)

  clust <- cluster::pam(sigMat, k = nClusters)

  GOAnnot <- clustToGO(clust, sigMat)

  gg <- plotGOHeatmap(sigMat, clust, GOAnnot)


  return(gg)
}

# Helper functions -------------------------------------------------------------

collapseSigGenes <- function(dgeList, p, genes = .data$Newid){
  genes <- rlang::enquo(genes)
  sigG <- purrr::map(dgeList, \(df) dplyr::filter(df, .data$PValue < p)) |>
    purrr::map(\(df) dplyr::select(df, !!genes)) |>
    dplyr::bind_rows() |>
    dplyr::distinct() |>
    unlist()

  return(sigG)
}

clustToGO <- function(clustering, sigMat, species = "Mouse", geneID = "symbol", maxTerms = 5){
  if(species != "Mouse"){
    stop("Currently, only Mouse is supported. Please submit an issue to add other species")
  }

  TableList <- clustering$clustering |>
    tibble::enframe(name = "Gene", value = "Cluster") |>
    dplyr::group_by(.data$Cluster) |>
    tidyr::nest() |>
    dplyr::mutate(geneUniverse = purrr::map(.data$data, \(data) factor(as.integer(rownames(sigMat) %in% data$Gene))),
           geneUniverse = purrr::map(.data$geneUniverse, \(univ) stats::setNames(univ, rownames(sigMat)))) |>
    dplyr::mutate(GoData = purrr::map(.data$geneUniverse, \(univ) methods::new("topGOdata",
                                                  ontology = "BP",
                                                  allGenes = univ,
                                                  annot = annFUN.org,
                                                  mapping = "org.Mm.eg.db",
                                                  ID = geneID))) |>
    dplyr::mutate(GoTest = purrr::map(.data$GoData, runTest, algorithm = "classic", statistic = "fisher")) |>
    dplyr::mutate(GoTable = purrr::map2(.data$GoData, .data$GoTest, \(data, test) GenTable(data, PValue = test, orderBy = "PValue", ranksOf = "PValue", topNodes = 20))) |>
    dplyr::select(.data$Cluster, .data$GoTable) |>
    dplyr::mutate(GoTable = purrr::map(.data$GoTable, \(df) dplyr::select(df, .data$GO.ID, .data$Term, .data$PValue))) |>
    tidyr::unnest(.data$GoTable) |>
    dplyr::group_by(.data$Cluster, .data$PValue) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::mutate(TermString = paste0(.data$GO.ID, ": ", .data$Term, " (P = ", .data$PValue, ")")) |>
    dplyr::group_by(.data$Cluster) |>
    dplyr::mutate(TermString = dplyr::case_when(min(.data$PValue) >= 0.05 ~ "No terms significantly enriched",
                                  TRUE ~ .data$TermString)) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$TermString == "No terms significantly enriched" | .data$PValue < 0.05) |>
    dplyr::select(.data$Cluster, .data$TermString) |>
    dplyr::distinct() |>
    dplyr::group_by(.data$Cluster) |>
    dplyr::slice_head(n = maxTerms) |>
    dplyr::summarise(TermString = paste(.data$TermString, collapse = "\n")) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = .data$Cluster, values_from = .data$TermString) |>
    unlist()

  return(TableList)
}

plotGOHeatmap <- function(matrix, clustering, GOAnnot){
  GOAnnot <- GOAnnot

  ht <- ComplexHeatmap::Heatmap(matrix,
                name = "zScore",
                cluster_columns = FALSE,
                show_row_dend = FALSE,
                show_row_names = FALSE,
                split = clustering$clustering,
                row_title = "@{GOAnnot[ x[1] ]}",
                row_title_rot = 0,
                gap = ggplot2::unit(2.5, "mm"),
                row_title_gp = grid::gpar(fontsize = 8))

  gg <- grid::grid.grabExpr(ComplexHeatmap::draw(ht)) |>
    ggpubr::as_ggplot()

  return(gg)
}
