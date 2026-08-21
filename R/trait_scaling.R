##########

#' Standard deviation of each trait across species
#'
#' Internal helper.  Replicate observations are averaged within a species first,
#' so that a heavily sampled species does not dominate the spread.
#'
#' @param traits A data frame whose first column is `species` and whose
#'   remaining columns are traits.
#' @return Named numeric vector, one entry per trait.
#' @noRd
trait_standard_deviations <- function(traits){

  trait_columns <- setdiff(colnames(traits), "species")

  values <- as.matrix(traits[, trait_columns, drop = FALSE])

  species <- as.character(traits$species)

  species_means <- vapply(X = seq_along(trait_columns),
                          FUN = function(column){

                            means <- tapply(X = values[, column],
                                            INDEX = species,
                                            FUN = mean,
                                            na.rm = TRUE)

                            stats::sd(as.numeric(means), na.rm = TRUE)

                          },
                          FUN.VALUE = numeric(1))

  stats::setNames(species_means, trait_columns)

}

#' Per-trait multipliers for the distance calculation
#'
#' Internal helper.  Branch lengths are Euclidean distances across traits, so
#' each trait enters in whatever units it was measured in and the trait with the
#' largest numbers dominates. These multipliers are applied to the reconstructed
#' values before the distance is taken, which is equivalent to rescaling the
#' trait data beforehand: Brownian motion reconstruction is equivariant under a
#' diagonal rescaling of the traits.
#'
#' @param traits A data frame whose first column is `species`.
#' @param scale_traits If TRUE, divide each trait by its standard deviation
#'   across species, so that traits contribute comparably regardless of units.
#' @param weights Optional weights on the squared differences, one per trait.
#' @return Numeric vector of multipliers, one per trait.
#' @noRd
trait_scaling_multipliers <- function(traits, scale_traits = FALSE, weights = NULL){

  trait_columns <- setdiff(colnames(traits), "species")

  multipliers <- rep(1, length(trait_columns))

  if(!is.logical(scale_traits) || length(scale_traits) != 1 || is.na(scale_traits)){

    stop("scale_traits must be TRUE or FALSE.", call. = FALSE)

  }

  if(scale_traits){

    deviations <- trait_standard_deviations(traits)

    usable <- is.finite(deviations) & deviations > 0

    if(any(!usable)){

      warning("Cannot scale ", sum(!usable), " of ", length(deviations),
              " traits, which have no variation across species; ",
              "they are left in their original units.",
              call. = FALSE)

    }

    multipliers[usable] <- 1/deviations[usable]

  }

  if(!is.null(weights)){

    if(!is.numeric(weights) || length(weights) != length(trait_columns) ||
       anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)){

      stop("weights must be a numeric vector of non-negative, finite values, ",
           "one per trait (", length(trait_columns), " here).",
           call. = FALSE)

    }

    if(all(weights == 0)){

      stop("weights must not be all zero.", call. = FALSE)

    }

    #Weights apply to the squared differences, so they enter as square roots.
    multipliers <- multipliers * sqrt(weights)

  }

  multipliers

}

#' Warn when one trait's units will dominate the distance
#'
#' Internal helper.  Only relevant when the traits are being used as supplied:
#' if the caller has asked for scaling or given weights, they have already made
#' the choice.
#'
#' @param traits A data frame whose first column is `species`.
#' @param ratio_threshold Warn when the largest and smallest trait standard
#'   deviations differ by more than this factor.
#' @return Invisibly NULL, called for the warning.
#' @noRd
warn_dominant_trait_scale <- function(traits, ratio_threshold = 10){

  deviations <- trait_standard_deviations(traits)

  usable <- deviations[is.finite(deviations) & deviations > 0]

  if(length(usable) < 2){

    return(invisible(NULL))

  }

  ratio <- max(usable)/min(usable)

  if(ratio <= ratio_threshold){

    return(invisible(NULL))

  }

  warning("Branch lengths are Euclidean distances across traits measured in their ",
          "original units, and these traits differ in spread by a factor of ",
          format(ratio, digits = 3), " (",
          names(which.max(usable)), " vs ", names(which.min(usable)),
          "), so the distance is dominated by the former. ",
          "Consider scale_traits = TRUE, or supplying weights.",
          call. = FALSE)

  invisible(NULL)

}

#' Apply the per-trait multipliers to a matrix of node values
#'
#' Internal helper.
#'
#' @param node_values Matrix of trait values, nodes by traits.
#' @param multipliers Numeric vector, one per trait.
#' @return The rescaled matrix.
#' @noRd
apply_trait_multipliers <- function(node_values, multipliers){

  if(all(multipliers == 1)){

    return(node_values)

  }

  sweep(x = node_values, MARGIN = 2, STATS = multipliers, FUN = "*")

}
