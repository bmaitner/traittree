##########

#' Check that PhyloMeasures is available
#'
#' Internal helper.  `PhyloMeasures` provides the phylogenetic diversity queries
#' and is only needed by the functions that compute them, so it is a suggested
#' rather than a required dependency: it was archived from CRAN in 2022 and has
#' to be installed from the archive, which should not stand between a user and
#' the rest of the package.
#'
#' The `installed` argument exists so that the failure message can be tested
#' without hiding the package.
#'
#' @param installed Whether the package is available. Defaults to checking.
#' @return Invisibly TRUE, or an error naming how to install it.
#' @noRd
require_phylomeasures <- function(installed = requireNamespace("PhyloMeasures", quietly = TRUE)){

  if(!installed){

    stop("PhyloMeasures is needed to compute the phylogenetic diversity metrics, ",
         "but it is not installed. It was archived from CRAN, so it has to come ",
         "from the archive:\n",
         "  remotes::install_url(\"https://cran.r-project.org/src/contrib/Archive/",
         "PhyloMeasures/PhyloMeasures_2.1.tar.gz\")\n",
         "It is compiled from source, so build tools are required (Rtools on ",
         "Windows, Xcode command line tools on macOS).",
         call. = FALSE)

  }

  invisible(TRUE)

}
