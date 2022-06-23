#' Convert DEG results tables to something Pathview can work with
#'
#' Helper function to take the output table from running dge_OneFactor(), filter
#' on an adjusted P value cutoff, and format it to something pathview can use.
#' Also can optionally enforce sentence case for mouse gene IDs (useful when you
#' get data from an informatics core that defaults to human-style gene names)
#' and deduplicating by selecting the most significant fold change when
#' duplicates are present.
#'
#' @param dgeRes Results table as returned by dge_OneFactor().
#' @param geneNames Tidy Eval selection of the column in dge_OneFactor that
#'   contains the gene names. Usually Newid if you converted names to MGI
#'   symbols and GeneID if you did not.
#' @param pval Numeric giving the P Value cutoff to use. Will keep only rows
#'   with adjusted P value *less than* this.
#' @param enforceSentenceCase Boolean; should gene names be forced into sentence
#'   case (capital first letter only)?
#' @param deduplicate Boolean; should genes be de-duplicated by keeping only the
#'   most significant instance of each one?
#'
#' @return A dataframe suitable for plotting KEGG pathways with plotKEGG().
#'
#' @importFrom rlang .data
#' @export
#'
dge_PrepKeggData <- function(dgeRes, geneNames = .data$Newid, pval = 0.05, enforceSentenceCase = FALSE, deduplicate = FALSE){
  if(!methods::is(dgeRes, "data.frame")){
    stop("dgeRes must be a single dataframe. Please vectorize if you want to use
         all outputs from dge_OneFactor, and make sure you don't include DGEList
         objects.")
  }

  df <- dplyr::filter(dgeRes, .data$PValue < pval)

  geneNames <- rlang::enquo(geneNames)

  if(enforceSentenceCase == TRUE){
    df <- dplyr::mutate(df, geneNames = stringr::str_to_sentence(!!geneNames))
  }

  if(deduplicate == TRUE){
    df <- df |>
      dplyr::group_by(!!geneNames) |>
      dplyr::filter(.data$PValue == min(.data$PValue)) |>
      dplyr::ungroup()
  }

  df <- dplyr::select(df, !!geneNames, .data$logFC) |>
    dplyr::filter(!is.na(!!geneNames)) |>
    tibble::column_to_rownames(var = rlang::quo_name(geneNames))

  return(df)
}

#' @examples
#' dgeRes <- system.file("extdata", "dgeRes_full.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'
#'   keggData <- dge_PrepKeggData(dgeRes, enforceSentenceCase = TRUE, deduplicate = TRUE)
