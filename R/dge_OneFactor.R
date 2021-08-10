#' Perform Differential Gene Expression analysis on a one-factor experiment
#'
#' This function handles a "default" call of edgeR to perform differential gene
#' expression testing on a count matrix file. The experiment must be a
#' one-factor design (this function isn't set up to handle multi-factor problems
#' yet), but can have any number of levels in that design, as long as they are
#' compared to a single reference group. Note that edgeR itself can handle much
#' more complex designs through limma, but having the ability to catch that
#' complexity in this function call is not so straightforward, and so very
#' complex experiments should probably be run more manually than this.
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
#' @param convertNames Boolean specifying whether names should be converted.
#'   Note that this requires an internet connection, so keep false if you are
#'   offline. Typically, this is used to convert the stable Ensembl ID that
#'   GeneWiz gives us into a more useful MGI symbol for downstream analysis.
#' @param convertFrom a string specifying a biomart object to convert the names
#'   from. Popular choices are "entrezgene_id" and "ensembl_gene_id".
#' @param convertTo a string specifying a biomart object (or character vector
#'   for several) to convert names to. Popular choices are "mgi_symbol".

#' @export
dge_OneFactor <- function(HitCountsMatrix,
                          Groups,
                          refGroup = NULL,
                          convertNames = FALSE,
                          convertFrom = NULL,
                          convertTo = NULL){

  if(length(Groups) != ncol(HitCountsMatrix)){
    stop("Length of Groups factor must match number of columns (samples) in
         HitCountsMatrix")
  }

  y <- edgeR::DGEList(counts = HitCountsMatrix, group = Groups)
  keep <- edgeR::filterByExpr(y)
  y <- y[keep,,keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  design <- stats::model.matrix(~ 0 + Groups)
  colnames(design) <- levels(Groups)

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design)

  if(is.null(refGroup)){
    refGroup <- levels(Groups[[1]])
  }

  compGroups <- levels(Groups)[levels(Groups) != refGroup]
  names(compGroups) <- as.character(compGroups)

  conlist <- purrr::map(compGroups,
                        ~limma::makeContrasts(contrasts = paste(..1,
                                                                "-",
                                                                refGroup,
                                                                sep = ""),
                                              levels = design))

  FCGeneUniverse <- purrr::map(conlist,
                               ~edgeR::glmTreat(fit, contrast = ..1)) |>
    purrr::map(~purrr::pluck(..1, "table")) |>
    purrr::map(~tibble::rownames_to_column(..1, var = "GeneID")) |>
    purrr::map(~dplyr::mutate(..1, GeneID = as.character(GeneID)))

  if(convertNames == TRUE){
    FCGeneUniverse <- purrr::map(FCGeneUniverse,
                                 ~convertBiomart(..1,
                                                 convert_from = convertFrom,
                                                 convert_to = convertTo))
  }

  if(length(FCGeneUniverse) == 1){
    FCGeneUniverse <- purrr::flatten(FCGeneUniverse)
  }

  return(FCGeneUniverse)
}




