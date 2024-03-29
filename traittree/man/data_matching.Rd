% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_matching_phylos_and_occurrences.R
\name{data_matching}
\alias{data_matching}
\title{Prune trait and phylogeny to a common set of species}
\usage{
data_matching(phylogeny, occurrences)
}
\arguments{
\item{phylogeny}{A phylogeny with time-scaled branch lengths}

\item{occurrences}{A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.}
}
\value{
List the pruned phylogeny and occurrences data
}
\description{
data_matching
}
\examples{
\dontrun{
Write example text
}
}
