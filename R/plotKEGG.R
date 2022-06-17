#' Create a KEGG pathway diagram with DEG data overlay.
#'
#' @param keggData A dataframe containing MGI symbol-formatted gene IDs,
#' @param pathway A string giving the pathway ID for the pathway of interest, in
#'   the format "00001". No species code is required.
#' @param species A string specifying the species to use. Currently, only
#'   "Mouse" is supported.
#' @param inputDir A path to a directory where pathview can download the KEGG
#'   files. Directory must exist before running this!
#' @param plotName An optional string to append to the output image file. Useful
#'   when vectorizing so each subsequent output doesn't overwrite the previous
#'   one.
#'
#' @return Nothing - instead, writes the output png to the working directory.
#'   Would love to have it not do that but this seems to be a feature of
#'   pathview.
#' @export
#'
plotKEGG <- function(keggData, pathway, species = "Mouse", inputDir = "KEGG Inputs/", plotName = NULL){
  if(species != "Mouse"){
    stop("plotKEGG currently only supports mouse data. Please submit an issue if you want additional species added.")
  }else{
    species <- "mmu"
  }
  pathview::pathview(gene.data = keggData,
           species = species,
           kegg.dir = inputDir,
           pathway.id = pathway,
           gene.idtype = "SYMBOL",
           same.layer = TRUE,
           multi.state = TRUE,
           out.suffix = plotName,
           low = list(gene = "cyan", cpd = "blue"),
           mid = list(gene = "white", cpd = "gray"),
           high = list(gene = "firebrick", cpd = "yellow"))
}

#' @examples
#' \dontrun{
#'   dgeRes <- system.file("extdata", "dgeRes_full.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#'   keggData <- dge_PrepKeggData(dgeRes, enforceSentenceCase = TRUE, deduplicate = TRUE)
#'
#'   plotKEGG(keggData, "04060", plotName = "Example")
#' }
