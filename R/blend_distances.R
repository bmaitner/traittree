##########

#' Coerce an input to a species-by-species distance matrix
#'
#' Internal helper.  Accepts a phylogeny, a `dist` object, an existing square
#' distance matrix, or a species-by-trait matrix, and returns a square matrix
#' with species names on both margins.
#'
#' @param x A phylogeny, `dist`, square distance matrix, or species-by-trait matrix.
#' @return A square numeric matrix with dimnames.
#' @noRd
as_distance_matrix <- function(x){

  if(inherits(x, "phylo")){ return(ape::cophenetic.phylo(x)) }

  if(inherits(x, "dist")){ return(as.matrix(x)) }

  x <- as.matrix(x)

  # a square, symmetric matrix is taken as distances already; anything else is
  # treated as species x traits and converted
  if(nrow(x) == ncol(x) &&
     isTRUE(all.equal(x, t(x), check.attributes = FALSE))){ return(x) }

  as.matrix(stats::dist(x))

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
#' Which is preferable depends on the question being asked.  A trait-scaled
#' tree is often a better phylogenetic component here than a dated one, so
#' passing the output of [scale_branches_multidimensional()] as `phylo_dist` is
#' worth trying.
#'
#' @param phylo_dist A phylogeny, or a phylogenetic distance matrix.
#' @param trait_dist A trait distance matrix, or a species-by-trait matrix, in
#'   which case Euclidean distances are taken.
#' @param a Weight on the phylogenetic component, between 0 (traits only) and
#'   1 (phylogeny only).  There is no default: see [tune_blend_weight()] for one
#'   way to choose it.
#' @param p The p-norm.  Defaults to 2, the convention in Cadotte et al. (2013).
#'   `p = 1` gives a straight linear mixture.
#' @return A square distance matrix over the species shared by the two inputs.
#' @references Cadotte, M.W., Albert, C.H. and Walker, S.C. (2013) The ecology
#'   of differences: assessing community assembly with trait and evolutionary
#'   distances. Ecology Letters 16:1234-1244.
#' @seealso [tune_blend_weight()], [scale_branches_multidimensional()]
#' @export
#' @importFrom ape cophenetic.phylo
#' @importFrom stats dist cor
#' @examples
#' tree <- ape::rtree(20)
#' traits <- matrix(stats::rnorm(40), nrow = 20,
#'                  dimnames = list(tree$tip.label, c("t1", "t2")))
#' blended <- blend_distances(phylo_dist = tree, trait_dist = traits, a = 0.25)
blend_distances <- function(phylo_dist,
                            trait_dist,
                            a,
                            p = 2){

  stopifnot(length(a) == 1, is.finite(a), a >= 0, a <= 1, length(p) == 1, p > 0)

  phylo_matrix <- as_distance_matrix(phylo_dist)
  trait_matrix <- as_distance_matrix(trait_dist)

  # species are matched by NAME: a blend built on positional alignment is wrong
  # whenever the trait table and the tree are ordered differently, which is the
  # usual case
  species <- intersect(row.names(phylo_matrix), row.names(trait_matrix))

  if(length(species) < 3){
    stop("fewer than 3 species shared between phylo_dist and trait_dist")
  }

  phylo_matrix <- phylo_matrix[species, species]
  trait_matrix <- trait_matrix[species, species]

  phylo_matrix <- phylo_matrix / max(phylo_matrix)
  trait_matrix <- trait_matrix / max(trait_matrix)

  blended <- if(p == 1){

    (a * phylo_matrix) + ((1 - a) * trait_matrix)

  }else{

    ((a * phylo_matrix^p) + ((1 - a) * trait_matrix^p))^(1 / p)

  }

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
#' The estimate is stable but tends to sit slightly below the weight that would
#' have been best, so it errs towards the traits and away from the phylogeny.
#'
#' Inspect `mean_score` before trusting the returned weight.  The criterion only
#' has something to work with when the held-out trait is predictable from the
#' remaining traits, the phylogeny, or both.  When it is predictable from
#' neither, the score curve is flat and the weight that comes back is arbitrary
#' rather than estimated.  A flat curve is itself informative: it says the
#' measured traits carry no signal about one another and the tree carries none
#' about them.
#'
#' @param phylo_dist A phylogeny, or a phylogenetic distance matrix.
#' @param traits A species-by-trait matrix or data frame, with species as row
#'   names.
#' @param a_grid Candidate weights to search.
#' @param p The p-norm, passed to [blend_distances()].
#' @param method The fold criterion.  `"mantel"`, the default, correlates the
#'   blended distances with the held-out trait's distance matrix across all
#'   pairs.  `"uniqueness"` instead correlates nearest-neighbour distinctiveness
#'   with the held-out trait's nearest-neighbour uniqueness; in one dimension
#'   that quantity is the gap to the nearest value in a sorted vector, which
#'   makes it noticeably noisier, and it is retained only for comparison.
#' @return A list with the chosen weight `a`, the grid searched, the mean score
#'   for each candidate weight, and the per-fold scores.
#' @seealso [blend_distances()]
#' @export
#' @examples
#' tree <- ape::rtree(30)
#' traits <- matrix(stats::rnorm(90), nrow = 30,
#'                  dimnames = list(tree$tip.label, c("t1", "t2", "t3")))
#' fit <- tune_blend_weight(phylo_dist = tree, traits = traits)
#' fit$a
tune_blend_weight <- function(phylo_dist,
                              traits,
                              a_grid = seq(0, 1, by = 0.05),
                              p = 2,
                              method = c("mantel", "uniqueness")){

  method <- match.arg(method)

  phylo_matrix <- as_distance_matrix(phylo_dist)
  trait_values <- as.matrix(traits)

  if(ncol(trait_values) < 2){
    stop("tune_blend_weight() needs at least 2 traits; with one, `a` has no effect")
  }

  species <- intersect(row.names(phylo_matrix), row.names(trait_values))

  if(length(species) < 3){
    stop("fewer than 3 species shared between phylo_dist and traits")
  }

  phylo_matrix <- phylo_matrix[species, species]
  phylo_matrix <- phylo_matrix / max(phylo_matrix)
  trait_values <- trait_values[species, , drop = FALSE]

  lower <- lower.tri(phylo_matrix)

  fold_scores <- vapply(X = seq_len(ncol(trait_values)),
                        FUN = function(held_out){

    remaining <- as.matrix(stats::dist(trait_values[, -held_out, drop = FALSE]))
    remaining <- remaining / max(remaining)

    held <- as.matrix(stats::dist(trait_values[, held_out, drop = FALSE]))
    held <- held / max(held)

    vapply(X = a_grid,
           FUN = function(a){

             blended <- blend_distances(phylo_dist = phylo_matrix,
                                        trait_dist = remaining,
                                        a = a,
                                        p = p)

             if(method == "mantel"){

               stats::cor(blended[lower], held[lower])

             }else{

               stats::cor(nearest_neighbour_distance(blended),
                          nearest_neighbour_distance(held))

             }

           },
           FUN.VALUE = numeric(1))

  },
  FUN.VALUE = numeric(length(a_grid)))

  mean_score <- rowMeans(fold_scores, na.rm = TRUE)

  list(a = a_grid[which.max(mean_score)],
       a_grid = a_grid,
       mean_score = mean_score,
       fold_scores = fold_scores,
       method = method,
       p = p)

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
