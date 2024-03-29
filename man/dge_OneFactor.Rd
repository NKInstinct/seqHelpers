% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dge_OneFactor.R
\name{dge_OneFactor}
\alias{dge_OneFactor}
\title{Perform Differential Gene Expression analysis on a one-factor experiment}
\usage{
dge_OneFactor(
  HitCountsMatrix,
  Groups,
  refGroup = NULL,
  comps = NULL,
  convertNames = FALSE,
  convertFrom = NULL,
  convertTo = NULL,
  retainDGEList = FALSE,
  ...
)
}
\arguments{
\item{HitCountsMatrix}{A matrix of hit counts, where each row is a gene and
each column is a sample. Note that this should be de-duplicated before
running through this function; i.e. each gene should be UNIQUE, and if they
aren't you should add the non-unique reads together with dplyr before you
start in on this!}

\item{Groups}{A factor showing what group each sample belongs to. The length
of this factor must match the number of columns in the HitCountsMatrix
exactly}

\item{refGroup}{The control group that each other group in the Groups factor
will be compared against. If not specified (refGroup = NULL), defaults to
the first level in the factor.}

\item{comps}{A string or character vector giving the comparisons you want to
do, in the form of "B-A". Passing a named character vector will give the
output the same names.}

\item{convertNames}{Boolean specifying whether names should be converted.
Note that this requires an internet connection, so keep false if you are
offline. Typically, this is used to convert the stable Ensembl ID that
GeneWiz gives us into a more useful MGI symbol for downstream analysis.}

\item{convertFrom}{a string specifying a biomart object to convert the names
from. Popular choices are "entrezgene_id" and "ensembl_gene_id".}

\item{convertTo}{a string specifying a biomart object (or character vector
for several) to convert names to. Popular choices are "mgi_symbol".}

\item{retainDGEList}{Boolean specifying whether the DGEList object built by
edgeR should be retained and appended to the output (useful for making
heatmaps later)}

\item{...}{Additional arguments to pass to edgeR::filterByExpr. Most likely
choices are "min.count", which specifies the minimum number of reads in at
least some samples, and "min.total.count", which specifies the minimum
number of reads in total, for a gene to be considered.}
}
\value{
A dataframe or list of dataframes showing all genes with their
fold-change and p value for the specified comparison.
}
\description{
This function handles a "default" call of edgeR to perform differential gene
expression testing on a count matrix file. The experiment must be a
one-factor design (this function isn't set up to handle multi-factor problems
yet), but can have any number of levels in that design. You can either
specify one reference group to compare all others against, or manually input
all comparisons you want it to do. See examples for details.
}
\examples{

  data <- system.file("extdata", "sampleMatrix.Rds",
                      package = "seqHelpers") |>
    readRDS()
  grouplist <- factor(c(rep("A", times = 3),
                        rep("B", times = 3),
                        rep("C", times = 3),
                        rep("D", times = 3)))

 # Default comparison, which will compare everything to "A"
 dge_OneFactor(data, Groups = grouplist)

 # Compare everything to D instead
 dge_OneFactor(data, Groups = grouplist, refGroup = "D")

 # Compare B to A and D to C only
 dge_OneFactor(data, Groups = grouplist, comps = c("B" = "B-A",
                                                   "D" = "D-C"))


}
