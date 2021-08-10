prepNormCounts <- function(counts, group){
  y <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(y)
  y <- y[keep,,keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  return(y)
}

prepFit <- function(y, group){
  design <- stats::model.matrix(~0+group)
  colnames(design) <- levels(group)

  fit <- edgeR::estimateDisp(y, design, robust = TRUE) |>
    edgeR::glmQLFit(design)

  prep <- list("fit" = fit,
               "design" = design)

  return(prep)
}

prepComps <- function(groups, ref, compare, data){
  if(is.null(ref) & is.null(compare)){
    ref <- levels(groups)[[1]]
  }

  if(is.null(compare)){
    compGroups <- levels(groups)[levels(groups) != ref]
    names(compGroups) <- as.character(compGroups)
    conlist <- purrr::map(compGroups,
                          ~limma::makeContrasts(contrasts = paste(..1,
                                                                  "-",
                                                                  ref,
                                                                  sep = ""),
                                                levels = data$design))
  } else {
    conlist <- purrr::map(compare,
                          ~limma::makeContrasts(contrasts = ..1,
                                                levels = data$design))
  }

  return(conlist)
}

doComps <- function(conlist, data){
  res <- purrr::map(conlist, ~edgeR::glmTreat(data$fit, contrast = ..1)) |>
    purrr::map(~purrr::pluck(..1, "table")) |>
    purrr::map(~tibble::rownames_to_column(..1, var = "GeneID")) |>
    purrr::map(~dplyr::mutate(..1, GeneID = as.character(GeneID)))

  return(res)
}
