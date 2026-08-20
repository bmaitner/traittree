##########

#' Identify species with replicate observations
#'
#' Internal helper.  Returns the names of species having two or more
#' non-missing values for at least one trait.
#'
#' Note that this counts observations per species *per trait*, not rows per
#' species.  Trait data are often stored with one row per species per source, so
#' a species can easily occupy several rows while still contributing only a
#' single observation of each trait; that species carries no information about
#' within-species variation.
#'
#' @param traits A data frame whose first column is `species` and whose
#'   remaining columns are traits.
#' @return Character vector of species names, possibly empty.
#' @noRd
replicated_species <- function(traits){

  trait_columns <- setdiff(colnames(traits), "species")

  if(nrow(traits) == 0 || length(trait_columns) == 0){

    return(character(0))

  }

  observed <- !is.na(traits[, trait_columns, drop = FALSE])

  observations_per_species <- rowsum(x = observed * 1,
                                     group = as.character(traits$species),
                                     reorder = FALSE)

  row.names(observations_per_species)[apply(observations_per_species > 1, MARGIN = 1, FUN = any)]

}

#' Does any species have replicate observations of a trait?
#'
#' Internal helper.
#'
#' @param traits A data frame whose first column is `species` and whose
#'   remaining columns are traits.
#' @return Single logical value.
#' @noRd
has_replicate_observations <- function(traits){

  length(replicated_species(traits)) > 0

}

#' Decide whether to estimate within-species (phenotypic) error
#'
#' Internal helper.
#'
#' `Rphylopars::phylopars()` has no declared default for `pheno_error`: leaving
#' the argument missing forces it to `TRUE` via the `pheno_correlated = TRUE`
#' default, whatever the data look like.  With a single observation per species
#' that setting is inert, as `phylopars()` then never fits the term at all and
#' returns no `phenocov`; reconstructions are identical either way.  Resolving
#' the value here therefore does not change results relative to the
#' `phylopars()` defaults.  What it does is state the specification explicitly
#' rather than inheriting it from an unrelated argument, and give callers a
#' documented way to override it.
#'
#' Where the setting does matter is with replicated data, for which the two
#' choices give substantially different reconstructions.  A single duplicated
#' row is enough to trip the switch, so sparse replication draws a warning.
#'
#' @param traits A data frame whose first column is `species` and whose
#'   remaining columns are traits.
#' @param pheno_error `NULL` to choose from the data, or `TRUE`/`FALSE` to force
#'   the setting.
#' @param sparse_replication_threshold Warn when replication is detected but the
#'   proportion of species carrying it falls below this value.
#' @return Single logical value, suitable for passing to `phylopars()`.
#' @noRd
resolve_pheno_error <- function(traits,
                                pheno_error = NULL,
                                sparse_replication_threshold = 0.05){

  if(!is.null(pheno_error)){

    if(!is.logical(pheno_error) || length(pheno_error) != 1 || is.na(pheno_error)){

      stop("pheno_error must be NULL, TRUE, or FALSE.", call. = FALSE)

    }

    return(pheno_error)

  }

  replicated <- replicated_species(traits)

  if(length(replicated) == 0){

    return(FALSE)

  }

  n_species <- length(unique(as.character(traits$species)))
  proportion_replicated <- length(replicated)/n_species

  if(proportion_replicated < sparse_replication_threshold){

    warning("Estimating within-species error from sparse replication: ",
            length(replicated), " of ", n_species, " species (",
            format(100 * proportion_replicated, digits = 2), "%) ",
            if(length(replicated) == 1) "has" else "have", " repeat observations. ",
            "A within-species error term will be fitted, which can change results ",
            "substantially and is considerably slower. If these repeats are duplicate ",
            "rows rather than genuine repeat measurements, pass pheno_error = FALSE.",
            call. = FALSE)

  }

  TRUE

}
