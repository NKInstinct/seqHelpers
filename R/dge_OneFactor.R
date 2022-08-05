#' Perform Differential Gene Expression analysis on a one-factor experiment
#'
#' This function handles a "default" call of edgeR to perform differential gene
#' expression testing on a count matrix file. The experiment must be a
#' one-factor design (this function isn't set up to handle multi-factor problems
#' yet), but can have any number of levels in that design. You can either
#' specify one reference group to compare all others against, or manually input
#' all comparisons you want it to do. See examples for details.
#'
#'
#' @param HitCountsMatrix A matrix of hit counts, where each row is a gene and
#'   each column is a sample. Note that this should be de-duplicated before
#'   running through this function; i.e. each gene should be UNIQUE, and if they
#'   aren't you should add the non-unique reads together with dplyr before you
#'   start in on this!
#' @param Groups A factor showing what group each sample belongs to. The length
#'   of this factor must match the number of columns in the HitCountsMatrix
#'   exactly
#' @param refGroup The control group that each other group in the Groups factor
#'   will be compared against. If not specified (refGroup = NULL), defaults to
#'   the first level in the factor.
#' @param comps A string or character vector giving the comparisons you want to
#'   do, in the form of "B-A". Passing a named character vector will give the
#'   output the same names.
#' @param convertNames Boolean specifying whether names should be converted.
#'   Note that this requires an internet connection, so keep false if you are
#'   offline. Typically, this is used to convert the stable Ensembl ID that
#'   GeneWiz gives us into a more useful MGI symbol for downstream analysis.
#' @param convertFrom a string specifying a biomart object to convert the names
#'   from. Popular choices are "entrezgene_id" and "ensembl_gene_id".
#' @param convertTo a string specifying a biomart object (or character vector
#'   for several) to convert names to. Popular choices are "mgi_symbol".
#' @param retainDGEList Boolean specifying whether the DGEList object built by
#'   edgeR should be retained and appended to the output (useful for making
#'   heatmaps later)
#' @param ... Additional arguments to pass to edgeR::filterByExpr. Most likely
#'   choices are "min.count", which specifies the minimum number of reads in at
#'   least some samples, and "min.total.count", which specifies the minimum
#'   number of reads in total, for a gene to be considered.
#'
#' @examples
#'
#'   data <- system.file("extdata", "sampleMatrix.Rds",
#'                       package = "seqHelpers") |>
#'     readRDS()
#'   grouplist <- factor(c(rep("A", times = 3),
#'                         rep("B", times = 3),
#'                         rep("C", times = 3),
#'                         rep("D", times = 3)))
#'
#'  # Default comparison, which will compare everything to "A"
#'  dge_OneFactor(data, Groups = grouplist)
#'
#'  # Compare everything to D instead
#'  dge_OneFactor(data, Groups = grouplist, refGroup = "D")
#'
#'  # Compare B to A and D to C only
#'  dge_OneFactor(data, Groups = grouplist, comps = c("B" = "B-A",
#'                                                    "D" = "D-C"))
#'
#'
#' @return A dataframe or list of dataframes showing all genes with their
#'     fold-change and p value for the specified comparison.

#' @export
dge_OneFactor <- function(HitCountsMatrix,
                          Groups,
                          refGroup = NULL,
                          comps = NULL,
                          convertNames = FALSE,
                          convertFrom = NULL,
                          convertTo = NULL,
                          retainDGEList = FALSE,
                          ...){

  # Input QC -------------------------------------------------------------------
  if(length(Groups) != ncol(HitCountsMatrix)){
    stop("Length of Groups factor must match number of columns (samples) in
         HitCountsMatrix")
  }

  # Prepare data & fit ---------------------------------------------------------
    # (see utils.R for these functions)
  prep <- prepNormCounts(HitCountsMatrix, Groups, ...) |>
    prepFit(group = Groups)

  # Perform comparisons --------------------------------------------------------
  Results <- prepComps(Groups, refGroup, comps, prep) |>
    doComps(prep)

  # Convert names --------------------------------------------------------------
  if(convertNames == TRUE){
    Results <- purrr::map(Results,
                          ~convertBiomart(..1,
                                          convert_from = convertFrom,
                                          convert_to = convertTo))
  }

  # Delist structure if only one comparison ------------------------------------
  if(length(Results) == 1 & retainDGEList == FALSE){
    Results <- purrr::flatten(Results) |>
      tibble::as_tibble()
  }

  # Keep the DGEList for heatmapping -------------------------------------------
  if(retainDGEList == TRUE){
    Results <- c(Results, "DGEList" = list(prepNormCounts(HitCountsMatrix, Groups)))
  }

  return(Results)
}
