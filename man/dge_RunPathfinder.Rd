% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dge_RunPathfinder.R
\name{dge_RunPathfinder}
\alias{dge_RunPathfinder}
\title{Run pathfindR analysis on output from dge_OneFactor}
\usage{
dge_RunPathfinder(
  dgeRes,
  geneNames = .data$Newid,
  species = "Mouse",
  reportFolder = "pathfindR_Results"
)
}
\arguments{
\item{dgeRes}{Dataframe; the output of running dge_OneFactor()}

\item{geneNames}{Column ID to choose the name of the column containing gene
names (tidy evaluation).}

\item{species}{Character. Currently only "Mouse" is supported, but others can
be added as needed.}

\item{reportFolder}{Path to folder where the pathfindR report will be saved.}
}
\value{
A dataframe showing all significantly altered KEGG pathways and the
constituent genes that are up- and down-regulated.
}
\description{
Essentially just a wrapper around run_pathfindR(), which unfortunately also
means it'll produce an HTML report and there's no way to stop it.
}
\examples{

dgeRes <- system.file("extdata", "dgeRes.Rds",
                      package = "seqHelpers") |>
    readRDS()
 \dontrun{
    dge_RunPathfinder(dgeRes)
 }

}
