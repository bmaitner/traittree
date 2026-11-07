##########

#' Phylogenetic diversity and its standardized index
#'
#' Faith's phylogenetic diversity (PD) of each of a set of communities, and the
#' phylogenetic diversity index (PDI) that standardizes it against the diversity
#' expected of the same number of species drawn at random. The null mean and
#' standard deviation are exact, not estimated by randomization, so the index is
#' fast and carries no sampling error.
#'
#' @param tree A phylogeny with branch lengths, or a \code{multiPhylo} (or list)
#'   of phylogenies on the same species. Its tips are the species pool.
#' @param communities The communities to measure, in either of two forms. A
#'   matrix or data frame with one row per community and one column per
#'   species, named by species, in which any non-zero entry counts as present;
#'   or a data frame of two columns, species name then community, with one row
#'   per occurrence, the "tidy" format used elsewhere in this package. A data
#'   frame is read as the second form when it has exactly two columns and the
#'   first holds names. Communities are returned in the order of the rows of a
#'   table, or the order in which they first appear among occurrences.
#' @param standardize Logical. Defaults to TRUE. If FALSE only PD is returned,
#'   which avoids computing the null distribution.
#' @param rarefy_to Optional richness to rarefy every community to, giving
#'   values that can be compared between communities of different richness (see
#'   Details). Either a whole number of at least 2, or \code{"min"} for the
#'   smallest richness among communities of two or more species. Defaults to
#'   NULL, meaning no rarefied values are computed.
#'
#' @details
#' PD is the total branch length of the minimal subtree connecting the species
#' of a community. The path from that subtree to the root of the phylogeny is
#' not included, so two sister species have a PD equal to the sum of their two
#' branches, however deep in the tree they sit. This is sometimes called
#' unrooted PD, and is the definition used by \code{PhyloMeasures::pd.query()}
#' and by \code{picante::pd(include.root = FALSE)}.
#'
#' PDI is
#' \deqn{PDI = \frac{PD - \mu}{\sigma}}{PDI = (PD - mean) / sd}
#' where the mean and standard deviation are those of PD across every set of
#' the same number of species drawn, uniformly and without replacement, from
#' the tips of \code{tree} (Sandel 2018). Negative values indicate species more
#' closely related than expected by chance. Both moments are computed exactly
#' from the formulas of Nipperess and Matsen (2013); see
#' \code{\link{phylo_diversity_moments}}.
#'
#' The species pool is every tip of the phylogeny as supplied. Tips that occur
#' in no community still count towards it, so prune the phylogeny first if the
#' pool should be narrower.
#'
#' A community of fewer than two species has no subtree, and its PD and PDI are
#' zero. A community holding every species in the phylogeny has a PDI of zero,
#' there being only one such community to compare it with.
#'
#' \strong{Comparing communities of different richness.} The standardization
#' guarantees that randomly assembled communities score zero whatever their
#' richness. It does not make other values comparable: under the same degree of
#' phylogenetic clustering, a richer community receives a more extreme index
#' (Sandel 2018). \code{rarefy_to} gives values free of that dependence, by
#' rarefying each community to a common richness. No random subsets are drawn.
#' The rarefied PD is the exact expected PD of a random subset of that many of
#' the community's species (Nipperess and Matsen 2013), and the rarefied index
#' is that expectation standardized against the null for the target richness,
#' which is the average of the index over every possible rarefaction.
#' Communities with fewer species than the target are NA. Rarefaction discards
#' information, more of it the lower the target, so the rarefied index
#' discriminates less well than the plain one; choose the largest target that
#' keeps the communities of interest.
#'
#' When several phylogenies are supplied and they share a topology, as the
#' trees returned by
#' \code{\link{scale_branches_multidimensional_with_variation}} do, they are
#' measured together at little more than the cost of one.
#'
#' @return For a single phylogeny, a data frame with one row per community and
#'   columns \code{community}, \code{richness} (the number of its species found
#'   in the phylogeny), \code{pd} and, unless \code{standardize = FALSE},
#'   \code{pdi}. With \code{rarefy_to}, also \code{pd_rarefied} and
#'   \code{pdi_rarefied}, and the richness used as the attribute
#'   \code{"rarefied_to"}.
#'
#'   For several phylogenies, a list of \code{community}, \code{richness}, and
#'   the same quantities as matrices with one row per community and one column
#'   per phylogeny, plus \code{rarefied_to} when used.
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#'   diversity. Biological Conservation 61: 1-10.
#'   \doi{10.1016/0006-3207(92)91201-3}
#'
#'   Nipperess, D. A. and Matsen, F. A. (2013). The mean and variance of
#'   phylogenetic diversity under rarefaction. Methods in Ecology and Evolution
#'   4: 566-572. \doi{10.1111/2041-210X.12042}
#'
#'   Sandel, B. (2018). Richness-dependence of phylogenetic diversity indices.
#'   Ecography 41: 837-844. \doi{10.1111/ecog.02967}
#'
#'   Tsirogiannis, C. and Sandel, B. (2016). PhyloMeasures: a package for
#'   computing phylogenetic biodiversity measures and their statistical moments.
#'   Ecography 39: 709-714. \doi{10.1111/ecog.01814}
#'
#' @seealso \code{\link{phylo_diversity_moments}} for the null mean and
#'   standard deviation themselves; \code{\link{calculate_pd_metric_rasters}}
#'   to map these quantities.
#'
#' @examples
#' # communities as occurrences: species, then the cell each occurs in
#' head(example_occurrences)
#'
#' diversity <- phylo_diversity(tree = example_tree,
#'                              communities = example_occurrences)
#'
#' head(diversity)
#'
#' # the same communities as a presence/absence matrix
#' presence <- table(example_occurrences$cell, example_occurrences$species)
#'
#' from_matrix <- phylo_diversity(tree = example_tree, communities = presence)
#'
#' all.equal(from_matrix$pd, diversity$pd)
#'
#' # values that can be compared between communities of different richness
#' rarefied <- phylo_diversity(tree = example_tree,
#'                             communities = example_occurrences,
#'                             rarefy_to = 4)
#'
#' head(rarefied[rarefied$richness >= 4, ])
#'
#' # several trees at once: uncertainty in a trait-scaled phylogeny carried
#' # through to the diversity of each community
#' trees <- scale_branches_multidimensional_with_variation(tree = example_tree,
#'                                                         traits = example_traits,
#'                                                         n_trees = 20)
#'
#' across_trees <- phylo_diversity(tree = trees, communities = example_occurrences)
#'
#' dim(across_trees$pdi)
#'
#' # spread of each community's index across the posterior sample
#' summary(apply(across_trees$pdi, 1, stats::sd))
#' @export
phylo_diversity <- function(tree, communities, standardize = TRUE, rarefy_to = NULL){

  trees <- pd_check_trees(tree)

  if(!is.logical(standardize) || length(standardize) != 1 || is.na(standardize)){

    stop("standardize must be TRUE or FALSE.", call. = FALSE)

  }

  if(!is.null(rarefy_to) && !identical(rarefy_to, "min") &&
     !(is.numeric(rarefy_to) && length(rarefy_to) == 1 && !is.na(rarefy_to) &&
       rarefy_to >= 2 && rarefy_to == round(rarefy_to))){

    stop("rarefy_to must be NULL, \"min\", or a single whole number of at least 2.",
         call. = FALSE)

  }

  input <- pd_read_communities(communities)

  n_tips <- length(trees[[1]]$tip.label)

  #Species that are not in the phylogeny cannot be placed on it. Say so, since
  #a mismatch in names is easier to make than to notice.
  known <- input$species %in% trees[[1]]$tip.label

  if(any(!known)){

    unknown <- unique(input$species[!known])

    warning(length(unknown), " species in communities ",
            if(length(unknown) == 1) "is not a tip" else "are not tips",
            " of the phylogeny and ",
            if(length(unknown) == 1) "was" else "were",
            " ignored, for example: ",
            paste(utils::head(unknown, 3), collapse = ", "), ".",
            call. = FALSE)

  }

  community <- input$community[known]
  species <- input$species[known]

  n_communities <- length(input$names)

  target <- NULL

  if(!is.null(rarefy_to)){

    if(identical(rarefy_to, "min")){

      distinct <- !duplicated(data.frame(community, species))

      per_community <- tabulate(community[distinct], nbins = n_communities)

      eligible <- per_community[per_community >= 2]

      target <- if(length(eligible) > 0) min(eligible) else 2

    }else{

      target <- as.integer(rarefy_to)

    }

    if(target > n_tips){

      stop("rarefy_to (", target, ") is larger than the number of species in the ",
           "phylogeny (", n_tips, ").", call. = FALSE)

    }

  }

  result <- pd_by_community(trees = trees,
                            community = community,
                            species = species,
                            n_communities = n_communities,
                            rarefy_to = target,
                            standardize = standardize)

  if(!is.null(target) && all(is.na(result$rarefied_pd))){

    warning("No community has ", target, " or more species, so the rarefied values ",
            "are NA throughout. Choose a smaller rarefy_to.", call. = FALSE)

  }

  #One phylogeny: a data frame, a row per community.
  if(inherits(tree, "phylo")){

    out <- data.frame(community = input$names,
                      richness = result$richness,
                      pd = result$pd[, 1],
                      stringsAsFactors = FALSE)

    if(standardize){

      out$pdi <- result$pdi[, 1]

    }

    if(!is.null(target)){

      out$pd_rarefied <- result$rarefied_pd[, 1]

      if(standardize){

        out$pdi_rarefied <- result$rarefied_pdi[, 1]

      }

      attr(out, "rarefied_to") <- target

    }

    return(out)

  }

  #Several: matrices, communities by phylogenies.
  label <- function(values){

    dimnames(values) <- list(input$names, names(trees))
    values

  }

  out <- list(community = input$names,
              richness = result$richness,
              pd = label(result$pd))

  if(standardize){

    out$pdi <- label(result$pdi)

  }

  if(!is.null(target)){

    out$pd_rarefied <- label(result$rarefied_pd)

    if(standardize){

      out$pdi_rarefied <- label(result$rarefied_pdi)

    }

    out$rarefied_to <- target

  }

  out

}

#' Null mean and standard deviation of phylogenetic diversity
#'
#' The exact mean and standard deviation of Faith's phylogenetic diversity (PD)
#' across every set of a given number of species drawn, uniformly and without
#' replacement, from the tips of a phylogeny. These are the quantities that
#' standardize PD into the phylogenetic diversity index, and the mean, as a
#' function of richness, is the rarefaction curve of PD.
#'
#' @param tree A phylogeny with branch lengths, or a \code{multiPhylo} (or list)
#'   of phylogenies on the same species.
#' @param richness The numbers of species to compute the moments for: whole
#'   numbers between 0 and the number of tips.
#'
#' @details
#' PD is the total branch length of the minimal subtree connecting a set of
#' species, without the path to the root; see \code{\link{phylo_diversity}}.
#'
#' With \eqn{s} tips in the phylogeny, \eqn{r} species drawn, \eqn{s_e} tips
#' below branch \eqn{e}, and \eqn{h(x) = {x \choose r} / {s \choose r}} the
#' probability that every species drawn falls within a given set of \eqn{x}
#' tips, a branch lies on the connecting subtree unless every species drawn is
#' below it or none is, so
#' \deqn{E[PD] = \sum_e w_e \, (1 - h(s_e) - h(s - s_e))}{E[PD] = sum_e w_e (1 - h(s_e) - h(s - s_e))}
#' where \eqn{w_e} is the length of the branch. The variance follows from the
#' probability that two branches both lie on the subtree, which takes one form
#' when one branch is an ancestor of the other and another when they are not
#' (Nipperess and Matsen 2013).
#'
#' Summed naively over pairs of branches the variance costs time quadratic in
#' the size of the phylogeny. Here branches are pooled by the number of tips
#' below them, and only ancestor-descendant pairs are treated individually,
#' following the approach of Tsirogiannis et al. (2014). The cost is of order
#' \eqn{n \log n} for a balanced phylogeny, rising towards quadratic for a
#' strongly unbalanced one, plus a term in the square of the number of distinct
#' clade sizes for each richness requested.
#'
#' None of these probabilities depends on branch lengths. Phylogenies that
#' share a topology are therefore handled together, at little more than the
#' cost of one.
#'
#' Fewer than two species have no subtree, so the mean and standard deviation
#' are zero. When every species is drawn there is one possible set, and the
#' standard deviation is exactly zero.
#'
#' @return For a single phylogeny, a data frame with columns \code{richness},
#'   \code{mean} and \code{sd}, one row per value of \code{richness}. For
#'   several, a list of \code{richness} and of \code{mean} and \code{sd} as
#'   matrices with one row per richness and one column per phylogeny.
#'
#' @references Nipperess, D. A. and Matsen, F. A. (2013). The mean and variance
#'   of phylogenetic diversity under rarefaction. Methods in Ecology and
#'   Evolution 4: 566-572. \doi{10.1111/2041-210X.12042}
#'
#'   Tsirogiannis, C., Sandel, B. and Kalvisa, A. (2014). New algorithms for
#'   computing phylogenetic biodiversity. Lecture Notes in Computer Science
#'   8701: 187-203. \doi{10.1007/978-3-662-44753-6_15}
#'
#' @seealso \code{\link{phylo_diversity}}
#'
#' @examples
#' # the rarefaction curve of PD for the example phylogeny, with its spread
#' moments <- phylo_diversity_moments(tree = example_tree, richness = 1:60)
#'
#' plot(moments$richness, moments$mean, type = "l",
#'      xlab = "species richness", ylab = "expected PD")
#' lines(moments$richness, moments$mean + 2 * moments$sd, lty = 2)
#' lines(moments$richness, moments$mean - 2 * moments$sd, lty = 2)
#'
#' # standardizing by hand reproduces the index
#' diversity <- phylo_diversity(example_tree, example_occurrences)
#'
#' null <- phylo_diversity_moments(example_tree, diversity$richness[1])
#'
#' (diversity$pd[1] - null$mean) / null$sd
#' diversity$pdi[1]
#' @export
phylo_diversity_moments <- function(tree, richness){

  trees <- pd_check_trees(tree)

  n_tips <- length(trees[[1]]$tip.label)

  if(!is.numeric(richness) || length(richness) == 0 || anyNA(richness) ||
     any(richness < 0) || any(richness > n_tips) || any(richness != round(richness))){

    stop("richness must be whole numbers between 0 and the number of tips (",
         n_tips, ").", call. = FALSE)

  }

  mean_pd <- matrix(data = NA_real_, nrow = length(richness), ncol = length(trees))
  sd_pd <- mean_pd

  #Phylogenies laid out like the first are handled together, any others singly.
  first <- trees[[1]]

  shared <- vapply(X = trees,
                   FUN = function(tree_i){

                     identical(tree_i$edge, first$edge) &&
                       identical(tree_i$tip.label, first$tip.label)

                   },
                   FUN.VALUE = logical(1))

  groups <- c(list(which(shared)), as.list(which(!shared)))

  for(group in groups){

    structure <- pd_tree_structure(trees[[group[1]]])

    lengths <- vapply(X = trees[group],
                      FUN = function(tree_i) as.numeric(tree_i$edge.length),
                      FUN.VALUE = numeric(structure$n_edges))

    sums <- pd_branch_sums(structure = structure,
                           lengths = matrix(data = lengths, nrow = structure$n_edges))

    for(value in unique(richness)){

      moments <- pd_null_moments(structure = structure, sums = sums, richness = value)

      rows <- which(richness == value)

      mean_pd[rows, group] <- rep(moments$mean, each = length(rows))
      sd_pd[rows, group] <- rep(moments$sd, each = length(rows))

    }

  }

  if(inherits(tree, "phylo")){

    return(data.frame(richness = richness,
                      mean = mean_pd[, 1],
                      sd = sd_pd[, 1]))

  }

  dimnames(mean_pd) <- list(richness, names(trees))
  dimnames(sd_pd) <- list(richness, names(trees))

  list(richness = richness, mean = mean_pd, sd = sd_pd)

}

#' Check and unpack the phylogenies given to the diversity functions
#'
#' Internal helper.
#'
#' @param tree A `phylo`, a `multiPhylo`, or a list of `phylo`.
#' @return A named list of `phylo` objects.
#' @noRd
pd_check_trees <- function(tree){

  if(inherits(tree, "phylo")){

    trees <- list(tree)

  }else if(inherits(tree, "multiPhylo")){

    #A multiPhylo may keep one copy of the tip labels in an attribute rather
    #than in each tree; put them back before treating it as a plain list.
    trees <- unclass(ape::.uncompressTipLabel(tree))

  }else if(is.list(tree)){

    trees <- tree

  }else{

    stop("tree must be a phylo object, or a multiPhylo or list of them.", call. = FALSE)

  }

  if(length(trees) == 0 || !all(vapply(trees, inherits, logical(1), "phylo"))){

    stop("tree must be a phylo object, or a multiPhylo or list of them.", call. = FALSE)

  }

  for(tree_i in trees){

    if(is.null(tree_i$edge.length) || anyNA(tree_i$edge.length)){

      stop("Every phylogeny needs a length for every branch.", call. = FALSE)

    }

    if(any(tree_i$edge.length < 0)){

      stop("Branch lengths must not be negative.", call. = FALSE)

    }

    if(anyDuplicated(tree_i$tip.label)){

      stop("Tip labels must be unique within a phylogeny.", call. = FALSE)

    }

  }

  #The species pool has to be the same for the phylogenies to be comparable.
  pool <- sort(trees[[1]]$tip.label)

  same_pool <- vapply(X = trees,
                      FUN = function(tree_i) identical(sort(tree_i$tip.label), pool),
                      FUN.VALUE = logical(1))

  if(!all(same_pool)){

    stop("Every phylogeny must have the same tips.", call. = FALSE)

  }

  if(is.null(names(trees))){

    names(trees) <- seq_along(trees)

  }

  trees

}

#' Read communities in either accepted form
#'
#' Internal helper.  Turns a communities-by-species table, or a two-column
#' table of occurrences, into one record per occurrence.
#'
#' @param communities See `phylo_diversity()`.
#' @return A list: `community`, an integer index for each occurrence;
#'   `species`, its species name; and `names`, the name of each community.
#' @noRd
pd_read_communities <- function(communities){

  #Two columns, names then community: one row per occurrence.
  tidy <- is.data.frame(communities) && ncol(communities) == 2 &&
    (is.character(communities[[1]]) || is.factor(communities[[1]]))

  if(tidy){

    if(anyNA(communities[[1]]) || anyNA(communities[[2]])){

      stop("communities contains missing values.", call. = FALSE)

    }

    species <- as.character(communities[[1]])
    key <- as.character(communities[[2]])

    community_names <- unique(communities[[2]])

    if(is.factor(community_names)){

      community_names <- as.character(community_names)

    }

    return(list(community = match(key, as.character(community_names)),
                species = species,
                names = community_names))

  }

  #Otherwise a table of communities by species.
  if(inherits(communities, "table") && length(dim(communities)) == 2){

    communities <- unclass(communities)

  }

  values <- as.matrix(communities)

  if(!(is.numeric(values) || is.logical(values))){

    stop("communities must be a numeric or logical table of communities by species, ",
         "or a two-column data frame of species and community.", call. = FALSE)

  }

  if(is.null(colnames(values))){

    stop("The columns of communities must be named by species.", call. = FALSE)

  }

  if(anyDuplicated(colnames(values))){

    stop("The columns of communities must have distinct species names.", call. = FALSE)

  }

  if(anyNA(values)){

    stop("communities contains missing values.", call. = FALSE)

  }

  if(any(values < 0)){

    stop("communities contains negative values.", call. = FALSE)

  }

  occupied <- which(values != 0, arr.ind = TRUE)

  community_names <- if(is.null(rownames(values))) seq_len(nrow(values)) else rownames(values)

  list(community = unname(occupied[, 1]),
       species = colnames(values)[occupied[, 2]],
       names = community_names)

}
