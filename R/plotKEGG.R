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
#' @param geneCols A 3-length list giving the low, medium, and high colour
#'   values to map to gene expression. Please keep them in that order, names are
#'   optional.
#' @param cpdCols As geneCols, but for compound values
#'
#' @return Nothing - instead, writes the output png to the working directory.
#'   Would love to have it not do that but this seems to be a feature of
#'   pathview.
#' @importFrom pathview pathview
#' @export
#'
plotKEGG <- function(keggData, pathway, species = "Mouse", inputDir = "KEGG Inputs/",
                     plotName = NULL, geneCols = list("low" = "cyan", "mid" = "white", "high" = "firebrick"),
                     cpdCols = list("low" = "blue", "mid" = "gray", "high" = "yellow")){
  if(species != "Mouse"){
    stop("plotKEGG currently only supports mouse data. Please submit an issue if you want additional species added.")
  }else{
    species <- "mmu"
  }
  pathview(gene.data = keggData,
           species = species,
           kegg.dir = inputDir,
           pathway.id = pathway,
           gene.idtype = "SYMBOL",
           same.layer = TRUE,
           multi.state = TRUE,
           out.suffix = plotName,
           low = list(gene = geneCols[[1]], cpd = cpdCols[[1]]),
           mid = list(gene = geneCols[[2]], cpd = cpdCols[[2]]),
           high = list(gene = geneCols[[3]], cpd = cpdCols[[3]]))
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
