##########

#' Warn when percent change is not a well behaved measure for these traits
#'
#' Internal helper.  Percent change is computed relative to the ancestral value
#' at the start of each branch.  That denominator is assumed to be strictly
#' positive and comfortably away from zero.  Neither holds automatically: values
#' near zero make the result explode, and negative values leave it without a
#' clear interpretation, since the ratio is then taken against a negative
#' baseline.  Both are routine for traits that have been log transformed or
#' centred, which is exactly when this measure is most likely to be reached for.
#'
#' The calculation is left alone; this only reports when its assumptions do not
#' hold, so the resulting branch lengths are not read at face value.
#'
#' @param denominator The ancestral values the change is divided by. A vector or
#'   a matrix of nodes by traits.
#' @param near_zero_fraction Values closer to zero than this fraction of a
#'   trait's own spread are counted as unstable.
#' @return Invisibly NULL, called for the warning.
#' @noRd
warn_unstable_percent_change <- function(denominator, near_zero_fraction = 0.01){

  denominator <- as.matrix(denominator)

  n_unstable <- 0L
  n_negative <- 0L

  for(column in seq_len(ncol(denominator))){

    values <- as.numeric(denominator[, column])

    spread <- stats::sd(values, na.rm = TRUE)

    if(is.finite(spread) && spread > 0){

      n_unstable <- n_unstable + sum(abs(values) < near_zero_fraction * spread, na.rm = TRUE)

    }

    n_negative <- n_negative + sum(values < 0, na.rm = TRUE)

  }

  if(n_unstable == 0 && n_negative == 0){

    return(invisible(NULL))

  }

  total <- length(denominator)

  details <- character(0)

  if(n_unstable > 0){

    details <- c(details,
                 paste0(n_unstable, " of ", total,
                        " are within ", 100 * near_zero_fraction,
                        "% of the trait's spread of zero, where the ratio is unstable"))

  }

  if(n_negative > 0){

    details <- c(details,
                 paste0(n_negative, " of ", total,
                        " are negative, where percent change has no clear interpretation"))

  }

  warning("percent = TRUE divides by the ancestral value at the start of each branch: ",
          paste(details, collapse = "; "),
          ". Consider percent = FALSE, or rescaling the traits to be strictly positive.",
          call. = FALSE)

  invisible(NULL)

}
