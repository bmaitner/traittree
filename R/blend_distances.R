##########

#' Coerce an input to a species-by-species distance matrix
#'
#' Internal helper.  Accepts a phylogeny, a `dist` object, an existing square
#' distance matrix, or a species-by-trait matrix, and returns a square matrix
#' with species names on both margins.
#'
#' @param x A phylogeny, `dist`, square distance matrix, or species-by-trait matrix.
#' @param argument The name of the user-facing argument `x` came from, used in
#'   error messages.
#' @return A square numeric matrix with dimnames.
#' @noRd
as_distance_matrix <- function(x, argument = "input"){

  if(inherits(x, "phylo")){ return(ape::cophenetic.phylo(x)) }

  if(inherits(x, "dist")){ return(as.matrix(x)) }

  x <- as.matrix(x)

  # a square input is read as distances only when it looks like distances on
  # every count: symmetric, zero on the diagonal and non-negative.  Squareness
  # and symmetry alone would silently misread a species-by-trait table that
  # happens to have as many traits as species.  Pass a `dist` object when the
  # distinction matters.
  looks_like_distances <- nrow(x) == ncol(x) &&
    isTRUE(all.equal(x, t(x), check.attributes = FALSE)) &&
    isTRUE(all(diag(x) == 0)) &&
    isTRUE(all(x >= 0))

  if(looks_like_distances){ return(name_both_margins(x, argument)) }

  if(is.null(row.names(x))){
    stop(argument, " needs species names as row names")
  }

  as.matrix(stats::dist(x))

}

##########

#' Give a square distance matrix names on both margins
#'
#' Internal helper.  Species are matched by name, so both margins have to carry
#' them; a distance matrix named on one margin only is named from the other.
#'
#' @param x A square numeric matrix.
#' @param argument The name of the user-facing argument `x` came from.
#' @return The matrix, with dimnames on both margins.
#' @noRd
name_both_margins <- function(x, argument){

  if(is.null(row.names(x)) && is.null(colnames(x))){
    stop(argument, " needs species names on at least one margin")
  }

  if(is.null(row.names(x))){ row.names(x) <- colnames(x) }

  if(is.null(colnames(x))){ colnames(x) <- row.names(x) }

  x

}

##########

#' Divide a distance matrix by its maximum
#'
#' Internal helper.  Scaling both sources to a maximum of one is what makes `a`
#' weight them against each other rather than against whatever units they were
#' measured in.  A matrix of all zeros (a constant trait, or a tree with no
#' branch lengths) has nothing to scale and is returned unchanged rather than
#' becoming `NaN`, which would otherwise poison even the endpoint blends.
#'
#' @param x A square numeric matrix.
#' @return The matrix, divided by its maximum when that maximum is positive.
#' @noRd
scale_to_unit_max <- function(x){

  largest <- max(x, na.rm = TRUE)

  if(!is.finite(largest) || largest <= 0){ return(x) }

  x / largest

}

##########

#' Check a blending weight
#'
#' Internal helper.
#'
#' @param a A candidate weight.
#' @return `TRUE`, invisibly; called for the error.
#' @noRd
check_blend_weight <- function(a){

  if(length(a) != 1 || !is.numeric(a) || !is.finite(a) || a < 0 || a > 1){
    stop("a must be a single number between 0 (traits only) and 1 (phylogeny only)")
  }

  invisible(TRUE)

}

##########

#' Check a grid of blending weights
#'
#' Internal helper.
#'
#' @param a_grid A candidate grid.
#' @return `TRUE`, invisibly; called for the error.
#' @noRd
check_blend_grid <- function(a_grid){

  if(!is.numeric(a_grid) || length(a_grid) < 1 || any(!is.finite(a_grid)) ||
     any(a_grid < 0) || any(a_grid > 1)){
    stop("a_grid must be a non-empty numeric vector of weights between 0 and 1")
  }

  invisible(TRUE)

}

##########

#' Check a p-norm
#'
#' Internal helper.  `p = Inf` is rejected rather than silently mishandled: in
#' that limit the blend is a component-wise maximum in which the weight `a`
#' drops out entirely, so it answers a different question.  `p < 1` is allowed
#' but warned about, since it breaks the triangle inequality.
#'
#' @param p A candidate p-norm.
#' @return `TRUE`, invisibly; called for the error and warning.
#' @noRd
check_p_norm <- function(p){

  if(length(p) != 1 || !is.numeric(p) || !is.finite(p) || p <= 0){
    stop("p must be a single finite positive number (p = Inf is not supported: ",
         "in that limit the blend is a component-wise maximum and a has no effect)")
  }

  if(p < 1){
    warning("p < 1 does not satisfy the triangle inequality, so the result is a ",
            "dissimilarity rather than a metric")
  }

  invisible(TRUE)

}

##########

#' Combine two prepared distance matrices
#'
#' Internal helper.  Assumes both matrices are already aligned and scaled, and
#' does no checking, so that the tuning loop can avoid re-coercing and
#' re-scaling its inputs at every point of the grid.
#'
#' @param phylo_matrix,trait_matrix Aligned, scaled square distance matrices.
#' @param a The weight on `phylo_matrix`.
#' @param p The p-norm.
#' @return A square distance matrix.
#' @noRd
blend_matrices <- function(phylo_matrix, trait_matrix, a, p){

  if(p == 1){

    (a * phylo_matrix) + ((1 - a) * trait_matrix)

  }else{

    ((a * phylo_matrix^p) + ((1 - a) * trait_matrix^p))^(1 / p)

  }

}

##########

#' Blend phylogenetic and trait distances
#'
#' Combines a phylogenetic distance matrix with a trait distance matrix into the
#' functional-phylogenetic distance of Cadotte, Albert and Walker (2013),
#'
#' \deqn{FPD = (a P^p + (1 - a) F^p)^{1/p}}
#'
#' where `P` and `F` are the phylogenetic and trait distances.  Each is divided
#' by its maximum first, so that `a` weights the two sources against each other
#' rather than against whatever units they happen to be measured in.
#'
#' This is the alternative to rescaling branch lengths: it combines the two
#' sources after the fact and returns a distance matrix, whereas
#' [scale_branches_multidimensional()] returns a phylogeny, on which
#' phylogenetic diversity and evolutionary distinctiveness remain defined.
#'
#' @section Choosing between the two routes:
#' The objective decides, more than the data do.  When the goal is to accumulate
#' diversity over a set of species (complementarity, surrogacy, reserve
#' selection), a trait-scaled tree did better than a blend in our simulations
#' and on real mammal and bird data.  When the goal is to rank individual
#' species by how distinctive they are (the EDGE and FUSE style of question),
#' the ordering reverses and a blend did better, up to about 75% trait coverage,
#' above which the measured traits alone are enough.  Blending also needs an
#' \eqn{O(n^2)} distance matrix and a weight, where rescaling is linear in tree
#' size and has no parameter; against that, a blend uses the measured trait
#' distances exactly, where rescaling passes them through an ancestral
#' reconstruction that shrinks them.
#'
#' A trait-scaled tree can be used as the phylogenetic component here, by
#' passing the output of [scale_branches_multidimensional()] as `phylo_dist`,
#' and that helped in simulation.  It did not help on real mammal data, where
#' the dated tree was the equal or better component at every trait count we
#' tried.  The value of a phylogenetic component in a blend appears to lie in
#' its independence from the measured traits, which rescaling gives up.  Worth
#' trying, in other words, but not worth assuming.
#'
#' @param phylo_dist A phylogeny, or a phylogenetic distance matrix.
#' @param trait_dist A trait distance matrix, or a species-by-trait matrix, in
#'   which case Euclidean distances are taken.  A square input is read as
#'   distances only if it is symmetric, zero on the diagonal and non-negative;
#'   pass a `dist` object to remove any ambiguity.
#' @param a Weight on the phylogenetic component, between 0 (traits only) and
#'   1 (phylogeny only).  There is no default: see [tune_blend_weight()] for one
#'   way to choose it.  On mammals and birds the best weight sat at 0.125 to
#'   0.25, with a flat score curve across that range, which is a reasonable
#'   place to start when the data cannot be cross-validated.
#' @param p The p-norm.  Defaults to 2, the convention in Cadotte et al. (2013).
#'   `p = 1` gives a straight linear mixture.  Fixing `p` rather than estimating
#'   it is deliberate: see [tune_blend_weight()].  Values below 1 are permitted
#'   but warn, since they break the triangle inequality and so return a
#'   dissimilarity rather than a metric.
#' @return A square distance matrix over the species shared by the two inputs.
#' @references Cadotte, M.W., Albert, C.H. and Walker, S.C. (2013) The ecology
#'   of differences: assessing community assembly with trait and evolutionary
#'   distances. Ecology Letters 16:1234-1244.
#' @seealso [tune_blend_weight()], [scale_branches_multidimensional()]
#' @export
#' @importFrom ape cophenetic.phylo
#' @importFrom stats dist cor sd
#' @examples
#' tree <- ape::rtree(20)
#' traits <- matrix(stats::rnorm(40), nrow = 20,
#'                  dimnames = list(tree$tip.label, c("t1", "t2")))
#' blended <- blend_distances(phylo_dist = tree, trait_dist = traits, a = 0.25)
blend_distances <- function(phylo_dist,
                            trait_dist,
                            a,
                            p = 2){

  check_blend_weight(a)
  check_p_norm(p)

  phylo_matrix <- as_distance_matrix(phylo_dist, "phylo_dist")
  trait_matrix <- as_distance_matrix(trait_dist, "trait_dist")

  # species are matched by NAME: a blend built on positional alignment is wrong
  # whenever the trait table and the tree are ordered differently, which is the
  # usual case
  species <- intersect(row.names(phylo_matrix), row.names(trait_matrix))

  if(length(species) < 3){
    stop("fewer than 3 species shared between phylo_dist and trait_dist")
  }

  phylo_matrix <- scale_to_unit_max(phylo_matrix[species, species])
  trait_matrix <- scale_to_unit_max(trait_matrix[species, species])

  blended <- blend_matrices(phylo_matrix = phylo_matrix,
                            trait_matrix = trait_matrix,
                            a = a,
                            p = p)

  dimnames(blended) <- list(species, species)

  blended

}

##########

#' Choose a blending weight by cross-validation
#'
#' The weight `a` in [blend_distances()] has no default and no established
#' estimator: the usual practice is to compute the blend across the whole range
#' of `a` and read the resulting curve as a diagnostic.  That is unhelpful when
#' a single value is needed.  This function chooses one from the data by holding
#' out each measured trait in turn and taking the weight whose blend best
#' reproduces the held-out trait's distance structure from the remaining traits
#' and the phylogeny.  Held-out measured traits stand in for the trait axes that
#' were never measured, which are what the blend is ultimately trying to
#' describe.
#'
#' At least two traits are required.  With one, holding it out leaves no trait
#' component and `a` has no effect on the result.
#'
#' @section How well this works:
#' Tested on mammals against the weight that would have been best for recovering
#' the full trait space, the returned weight landed within one grid step of that
#' weight in every subset tried and matched it exactly in about half, at a cost
#' of 0.04 in correlation.  It is biased low, selecting `a = 0` in roughly 85% of
#' subsets where the best weight was usually 0.125, so the penalty for taking the
#' argmax at face value is small only because the score curve is flat.  Read
#' `a_plateau` alongside `a`: it reports the range of weights the folds cannot
#' distinguish from the best one, and it is usually wide.
#'
#' Inspect `mean_score` before trusting the returned weight.  The criterion only
#' has something to work with when the held-out trait is predictable from the
#' remaining traits, the phylogeny, or both.  When it is predictable from
#' neither, the score curve is flat and the weight that comes back is arbitrary
#' rather than estimated.  A flat curve is itself informative: it says the
#' measured traits carry no signal about one another and the tree carries none
#' about them.
#'
#' @section Why p is fixed rather than tuned:
#' The argument that motivates estimating `a` would seem to apply to `p` as
#' well, which is simply set to the convention.  It does not hold up.  Scoring a
#' two-dimensional grid over both showed the achievable performance changing by
#' only 0.04 to 0.07 in correlation across `p` from 0.1 to 5, cross-validation
#' recovering the best `p` in 9% to 55% of subsets depending on the grid, and
#' tuning `p` alongside `a` beating a fixed `p = 2` in about a fifth of subsets
#' (0.806 against 0.794 on average).  `p` is weakly identified, so it is fixed
#' and `a` alone is tuned.
#'
#' @param phylo_dist A phylogeny, or a phylogenetic distance matrix.
#' @param traits A species-by-trait matrix or data frame, with species as row
#'   names.
#' @param a_grid Candidate weights to search: a non-empty numeric vector of
#'   values between 0 and 1.
#' @param p The p-norm, passed to [blend_distances()].
#' @param method The fold criterion.  `"mantel"`, the default, correlates the
#'   blended distances with the held-out trait's distance matrix across all
#'   pairs.  `"uniqueness"` instead correlates nearest-neighbour distinctiveness
#'   with the held-out trait's nearest-neighbour uniqueness; in one dimension
#'   that quantity is the gap to the nearest value in a sorted vector, which
#'   makes it noticeably noisier, and it is retained only for comparison.  On
#'   mammals the Mantel criterion matched the best weight in 55% of subsets
#'   against 36% for the other, at a lower cost (0.04 against 0.06).
#' @return A list with the chosen weight `a`, the range `a_plateau` of weights
#'   scoring within one standard error of it, the grid searched, the mean score
#'   and its standard error for each candidate weight, and the per-fold scores.
#' @seealso [blend_distances()]
#' @export
#' @examples
#' tree <- ape::rtree(30)
#' traits <- matrix(stats::rnorm(90), nrow = 30,
#'                  dimnames = list(tree$tip.label, c("t1", "t2", "t3")))
#' fit <- tune_blend_weight(phylo_dist = tree, traits = traits)
#' fit$a
#' fit$a_plateau
tune_blend_weight <- function(phylo_dist,
                              traits,
                              a_grid = seq(0, 1, by = 0.05),
                              p = 2,
                              method = c("mantel", "uniqueness")){

  method <- match.arg(method)

  check_blend_grid(a_grid)
  check_p_norm(p)

  phylo_matrix <- as_distance_matrix(phylo_dist, "phylo_dist")
  trait_values <- as.matrix(traits)

  if(ncol(trait_values) < 2){
    stop("tune_blend_weight() needs at least 2 traits; with one, `a` has no effect")
  }

  if(is.null(row.names(trait_values))){
    stop("traits needs species names as row names")
  }

  species <- intersect(row.names(phylo_matrix), row.names(trait_values))

  if(length(species) < 3){
    stop("fewer than 3 species shared between phylo_dist and traits")
  }

  phylo_matrix <- scale_to_unit_max(phylo_matrix[species, species])
  trait_values <- trait_values[species, , drop = FALSE]

  lower <- lower.tri(phylo_matrix)

  fold_scores <- vapply(X = seq_len(ncol(trait_values)),
                        FUN = function(held_out){

    remaining <- as.matrix(stats::dist(trait_values[, -held_out, drop = FALSE]))
    remaining <- scale_to_unit_max(remaining)

    held <- as.matrix(stats::dist(trait_values[, held_out, drop = FALSE]))
    held <- scale_to_unit_max(held)

    vapply(X = a_grid,
           FUN = function(a){

             blended <- blend_matrices(phylo_matrix = phylo_matrix,
                                       trait_matrix = remaining,
                                       a = a,
                                       p = p)

             if(method == "mantel"){

               scoring_correlation(blended[lower], held[lower])

             }else{

               scoring_correlation(nearest_neighbour_distance(blended),
                                   nearest_neighbour_distance(held))

             }

           },
           FUN.VALUE = numeric(1))

  },
  FUN.VALUE = numeric(length(a_grid)))

  # vapply drops to a vector when the grid holds a single candidate, which
  # rowMeans() cannot take
  fold_scores <- matrix(data = fold_scores,
                        nrow = length(a_grid),
                        dimnames = list(NULL, colnames(trait_values)))

  scored_folds <- rowSums(is.finite(fold_scores))
  mean_score <- rowMeans(fold_scores, na.rm = TRUE)
  mean_score[scored_folds == 0] <- NA_real_

  standard_error <- apply(X = fold_scores,
                          MARGIN = 1,
                          FUN = function(scores){
                            stats::sd(scores, na.rm = TRUE) /
                              sqrt(sum(is.finite(scores)))
                          })

  best <- which.max(mean_score)

  if(length(best) == 0){

    warning("no candidate weight could be scored: the held-out traits have no ",
            "variation, or nothing in the blend varies either")

    chosen <- NA_real_
    plateau <- c(NA_real_, NA_real_)

  }else{

    chosen <- a_grid[best]

    # the curve is usually flat, so report the weights the folds cannot tell
    # apart from the best one rather than the argmax alone
    tolerance <- standard_error[best]
    if(!is.finite(tolerance)){ tolerance <- 0 }

    indistinguishable <- which(is.finite(mean_score) &
                                 mean_score >= mean_score[best] - tolerance)

    plateau <- range(a_grid[indistinguishable])

  }

  list(a = chosen,
       a_plateau = plateau,
       a_grid = a_grid,
       mean_score = mean_score,
       standard_error = standard_error,
       fold_scores = fold_scores,
       method = method,
       p = p)

}

##########

#' Correlate two fold scores
#'
#' Internal helper.  A held-out trait that does not vary, or a blend that could
#' not be formed, has no correlation to report; `stats::cor()` warns and returns
#' `NA` in that case, and this returns the `NA` without the warning, since the
#' fold is dropped from the mean either way.
#'
#' @param x,y Numeric vectors.
#' @return A correlation, or `NA_real_`.
#' @noRd
scoring_correlation <- function(x, y){

  spreads <- c(stats::sd(x, na.rm = TRUE), stats::sd(y, na.rm = TRUE))

  if(any(!is.finite(spreads)) || any(spreads == 0)){ return(NA_real_) }

  stats::cor(x, y)

}

##########

#' Distance from each species to its nearest neighbour
#'
#' Internal helper.
#'
#' @param distance_matrix A square distance matrix.
#' @return Named numeric vector, one element per species.
#' @noRd
nearest_neighbour_distance <- function(distance_matrix){

  diag(distance_matrix) <- Inf

  apply(X = distance_matrix,
        MARGIN = 1,
        FUN = min)

}
