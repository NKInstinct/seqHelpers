#' @importFrom rlang .data
#'
convertBiomart <- function(data, convert_from = "ensemble_id", convert_to = "mgi_symbol"){
  ensembl <- biomaRt::useEnsembl("genes", dataset = "mmusculus_gene_ensembl")
  geneList <- data$GeneID
  geneTable <- biomaRt::getBM(attributes = c(convert_from, convert_to), values = geneList, mart = ensembl) |>
    dplyr::rename(GeneID = convert_from, Newid = dplyr::all_of(convert_to)) |>
    dplyr::mutate(GeneID = as.character(.data$GeneID))

  data <- dplyr::left_join(data, geneTable, by = "GeneID")
}
