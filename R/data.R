##########

#' Simulated phylogeny for the package examples
#'
#' A time-scaled phylogeny of 60 imaginary species, used throughout the examples
#' and the vignette. Branch lengths are scaled so that the root is at depth 1.
#'
#' The data are simulated, not observed. They exist so that the examples run
#' quickly and reproducibly, and should not be used for anything else. The code
#' that produced them is in \code{data-raw/generate_example_data.R}.
#'
#' @format An object of class \code{phylo} with 60 tips.
#' @seealso [example_traits], [example_occurrences]
"example_tree"

#' Simulated trait data for the package examples
#'
#' Three correlated traits for the 60 species in [example_tree], evolved under
#' Brownian motion and then placed on plausible measurement scales. The traits
#' are deliberately given different units, so that the effect of
#' \code{scale_traits} is visible: their standard deviations differ by roughly a
#' factor of eleven.
#'
#' The data are simulated, not observed. There are no missing values, though the
#' package's functions do accept them.
#'
#' @format A data frame with 60 rows and 4 columns:
#' \describe{
#'   \item{species}{Species name, matching the tip labels of [example_tree].}
#'   \item{body_mass_g}{Body mass, in grams.}
#'   \item{wing_mm}{Wing length, in millimetres.}
#'   \item{bill_mm}{Bill depth, in millimetres.}
#' }
#' @seealso [example_tree], [example_occurrences]
"example_traits"

#' Simulated occurrence data for the package examples
#'
#' Occurrences of the 60 species in [example_tree] across a 10 by 10 grid, in the
#' "tidy" format the package's spatial functions expect: one row per species per
#' occupied cell. Species occupy clustered sets of cells rather than being
#' scattered uniformly.
#'
#' The cell numbers correspond to a raster of 10 rows by 10 columns, which the
#' examples build with \code{raster::raster(nrows = 10, ncols = 10)}.
#'
#' The data are simulated, not observed.
#'
#' @format A data frame with 375 rows and 2 columns:
#' \describe{
#'   \item{species}{Species name, matching the tip labels of [example_tree].}
#'   \item{cell}{Raster cell number the species occurs in.}
#' }
#' @seealso [example_tree], [example_traits]
"example_occurrences"
