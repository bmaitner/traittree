##########

# Phylogenetic diversity (Faith's PD) and its standardized index (PDI), computed
# without external dependencies.
#
# PD here is the total branch length of the minimal subtree connecting a set of
# tips; the path from that subtree to the root is not included. PDI standardizes
# PD against every set of the same richness drawn uniformly, without
# replacement, from the tips of the tree (Sandel 2018, eq. 1):
#
#     PDI = (PD - mean) / sd
#
# The mean and variance are exact, from Nipperess and Matsen (2013, eq. 5 and
# Appendix). With s tips in the tree, r in the sample, s_e tips below edge e, and
# h(x) = choose(x, r) / choose(s, r) the probability that every draw falls in a
# given set of x tips:
#
#     P(e in the subtree)      = 1 - h(s_e) - h(s - s_e)
#     P(e and f), disjoint     = 1 - h(s - s_e) - h(s - s_f) + h(s - s_e - s_f)
#     P(e and f), e inside f   = 1 - h(s_f) - h(s - s_e) + h(s_f - s_e)
#
# and E[PD] and E[PD^2] are sums of branch lengths weighted by these.
#
# Summing over every pair of edges is quadratic in the size of the tree. Two
# things avoid that. The disjoint probability depends on the edges only through
# the sizes of their clades, so edges are pooled by clade size; and treating
# every pair as disjoint is wrong only for ancestor-descendant pairs, which are
# enumerated and corrected. There are as many of those as the summed depth of
# the edges, which is of order n log n for a balanced tree.
#
# None of the probabilities involves a branch length. Trees that share a
# topology, as every tree drawn by scale_branches_multidimensional_with_variation()
# does, therefore share all of the work except one weighted sum each.
#
# References:
#   Nipperess, D. A. and Matsen, F. A. (2013). The mean and variance of
#     phylogenetic diversity under rarefaction. Methods in Ecology and Evolution
#     4: 566-572.
#   Sandel, B. (2018). Richness-dependence of phylogenetic diversity indices.
#     Ecography 41: 837-844.
#   Tsirogiannis, C., Sandel, B. and Kalvisa, A. (2014). New algorithms for
#     computing phylogenetic biodiversity. LNCS 8701: 187-203.

#' Topology of a tree, arranged for phylogenetic diversity calculations
#'
#' Internal helper.  Everything here depends on the topology alone, so it is
#' computed once and reused for every set of branch lengths on that topology.
#'
#' @param tree A phylogeny.
#' @return A list: the number of tips and edges; for each tip, the edges on its
#'   path to the root; the number of tips below each edge and the distinct values
#'   of that number ("classes"); and every pair of edges in which one is an
#'   ancestor of the other, pooled by the clade sizes of the two.
#' @noRd
pd_tree_structure <- function(tree){

  edge <- tree$edge
  n_tips <- length(tree$tip.label)
  n_edges <- nrow(edge)

  #The edge above each edge, NA for edges that leave the root.
  parent_edge <- match(edge[, 1], edge[, 2])

  #Every (descendant, ancestor) pair of edges, found by walking all edges up the
  #tree together, one step per pass.
  descendants <- list()
  ancestors <- list()

  current <- seq_len(n_edges)
  above <- parent_edge

  repeat{

    keep <- !is.na(above)

    if(!any(keep)){

      break

    }

    current <- current[keep]
    above <- above[keep]

    descendants[[length(descendants) + 1]] <- current
    ancestors[[length(ancestors) + 1]] <- above

    above <- parent_edge[above]

  }

  descendants <- unlist(descendants)
  ancestors <- unlist(ancestors)

  if(is.null(descendants)){

    descendants <- integer(0)
    ancestors <- integer(0)

  }

  #The edges between each tip and the root: its own edge and that edge's
  #ancestors.
  tip_edge <- match(seq_len(n_tips), edge[, 2])

  from_tip <- edge[descendants, 2] <= n_tips

  path_edges <- split(x = c(tip_edge, ancestors[from_tip]),
                      f = factor(c(seq_len(n_tips), edge[descendants[from_tip], 2]),
                                 levels = seq_len(n_tips)))

  path_edges <- unname(path_edges)

  clade_size <- tabulate(unlist(path_edges), nbins = n_edges)

  #Pool edges by clade size, and ancestor-descendant pairs by the two sizes.
  class_sizes <- sort(unique(clade_size))
  edge_class <- match(clade_size, class_sizes)

  pair_key <- (edge_class[descendants] - 1) * length(class_sizes) + edge_class[ancestors]
  pair_keys <- sort(unique(pair_key))

  list(n_tips = n_tips,
       n_edges = n_edges,
       path_edges = path_edges,
       path_length = lengths(path_edges),
       clade_size = clade_size,
       class_sizes = class_sizes,
       edge_class = edge_class,
       descendants = descendants,
       ancestors = ancestors,
       pair_group = match(pair_key, pair_keys),
       group_lower = class_sizes[(pair_keys - 1) %/% length(class_sizes) + 1],
       group_upper = class_sizes[(pair_keys - 1) %% length(class_sizes) + 1])

}

#' Branch-length sums needed for the moments of phylogenetic diversity
#'
#' Internal helper.  Reduces one or more sets of branch lengths on a shared
#' topology to the handful of pooled sums the mean and variance are built from.
#'
#' @param structure Output of `pd_tree_structure()`.
#' @param lengths Matrix of branch lengths, edges by trees.
#' @return A list of the total length, the length and squared length pooled by
#'   clade size, and the products of ancestor-descendant pairs pooled by the two
#'   clade sizes. Each has one column per tree.
#' @noRd
pd_branch_sums <- function(structure, lengths){

  n_groups <- length(structure$group_lower)

  nested <- matrix(data = 0, nrow = n_groups, ncol = ncol(lengths))

  #Products of ancestor-descendant pairs, accumulated in pieces so that the
  #pairs-by-trees table never has to exist in full.
  n_pairs <- length(structure$descendants)

  if(n_pairs > 0){

    chunk <- max(1, floor(4e6/ncol(lengths)))

    for(start in seq(from = 1, to = n_pairs, by = chunk)){

      rows <- start:min(start + chunk - 1, n_pairs)

      pooled <- rowsum(x = lengths[structure$descendants[rows], , drop = FALSE] *
                         lengths[structure$ancestors[rows], , drop = FALSE],
                       group = structure$pair_group[rows],
                       reorder = FALSE)

      target <- as.integer(row.names(pooled))

      nested[target, ] <- nested[target, ] + pooled

    }

  }

  list(total = colSums(lengths),
       by_class = unname(rowsum(lengths, structure$edge_class)),
       squared_by_class = unname(rowsum(lengths^2, structure$edge_class)),
       nested = nested)

}

#' Mean and standard deviation of phylogenetic diversity for a given richness
#'
#' Internal helper.  The moments of PD over every set of `richness` tips drawn
#' uniformly without replacement from the tips of the tree.
#'
#' @param structure Output of `pd_tree_structure()`.
#' @param sums Output of `pd_branch_sums()`.
#' @param richness Single number of tips in the sample.
#' @return A list of `mean` and `sd`, each with one element per tree.
#' @noRd
pd_null_moments <- function(structure, sums, richness){

  n_trees <- length(sums$total)
  n_tips <- structure$n_tips

  #With fewer than two tips there is no subtree to measure.
  if(richness < 2){

    return(list(mean = rep(0, n_trees), sd = rep(0, n_trees)))

  }

  #h(x): the probability that every draw falls inside a given set of x tips.
  log_all <- lchoose(n_tips, richness)

  all_inside <- function(x){

    out <- numeric(length(x))
    possible <- x >= richness
    out[possible] <- exp(lchoose(x[possible], richness) - log_all)
    out

  }

  sizes <- structure$class_sizes

  inside <- all_inside(sizes)
  outside <- all_inside(n_tips - sizes)

  mean_pd <- as.vector((1 - inside - outside) %*% sums$by_class)

  #Every sample holds every tip, so there is nothing to vary.
  if(richness >= n_tips){

    return(list(mean = mean_pd, sd = rep(0, n_trees)))

  }

  #E[PD^2], first as though every pair of edges had disjoint clades...
  neither <- matrix(data = all_inside(n_tips - outer(sizes, sizes, "+")),
                    nrow = length(sizes))

  second_moment <- sums$total^2 -
    2 * sums$total * as.vector(outside %*% sums$by_class) +
    colSums(sums$by_class * (neither %*% sums$by_class))

  #...then corrected for the pairs where that is false: an edge with itself,
  #and an edge with each of its ancestors. The second kind appears twice in the
  #double sum.
  correction <- function(lower, upper){

    all_inside(n_tips - upper) - all_inside(upper) +
      all_inside(upper - lower) - all_inside(n_tips - lower - upper)

  }

  second_moment <- second_moment +
    as.vector(correction(sizes, sizes) %*% sums$squared_by_class)

  if(length(structure$group_lower) > 0){

    second_moment <- second_moment +
      2 * as.vector(correction(structure$group_lower, structure$group_upper) %*% sums$nested)

  }

  variance <- second_moment - mean_pd^2

  #The subtraction leaves rounding noise where the true variance is zero.
  variance[variance < 1e-14 * mean_pd^2] <- 0

  list(mean = mean_pd, sd = sqrt(variance))

}

#' Phylogenetic diversity of many communities at once
#'
#' Internal helper.  For each community, the total length of the edges that have
#' some but not all of its tips below them, which is the minimal subtree
#' connecting them.
#'
#' With `rarefy_to`, also the expected PD of a random subset of that many of the
#' community's tips, from the same formula as the null mean applied to the
#' community's own subtree.
#'
#' @param structure Output of `pd_tree_structure()`.
#' @param lengths Matrix of branch lengths, edges by trees.
#' @param community Integer community index of each occurrence.
#' @param tip Integer tip index of each occurrence. Each community-tip pair must
#'   appear once.
#' @param n_communities Number of communities.
#' @param rarefy_to Optional richness to rarefy to.
#' @param max_communities,max_work Limits on the size of a block: the number of
#'   communities, and the summed path lengths of their tips. The defaults keep
#'   memory use modest; they are arguments so that the blocking can be tested.
#' @return A list: `richness`, a vector; `pd`, a communities by trees matrix;
#'   and, if requested, `rarefied`, the same shape, NA where a community has
#'   fewer tips than `rarefy_to`.
#' @noRd
pd_of_communities <- function(structure, lengths, community, tip, n_communities,
                              rarefy_to = NULL,
                              max_communities = NULL,
                              max_work = NULL){

  n_edges <- structure$n_edges
  n_trees <- ncol(lengths)

  if(is.null(max_communities)){

    max_communities <- max(1, floor(4e6/n_edges))

  }

  if(is.null(max_work)){

    max_work <- 2e7/n_trees

  }

  in_order <- order(community)
  community <- community[in_order]
  tip <- tip[in_order]

  richness <- tabulate(community, nbins = n_communities)

  pd <- matrix(data = 0, nrow = n_communities, ncol = n_trees)

  rarefied <- NULL

  if(!is.null(rarefy_to)){

    rarefied <- matrix(data = NA_real_, nrow = n_communities, ncol = n_trees)

  }

  if(length(community) == 0){

    return(list(richness = richness, pd = pd, rarefied = rarefied))

  }

  last_row <- cumsum(richness)
  first_row <- last_row - richness + 1

  #Work in blocks of communities, sized so that neither the communities-by-edges
  #table of counts nor the table of included edges by trees grows large.
  work <- numeric(n_communities)
  per_community <- rowsum(structure$path_length[tip], community)
  work[as.integer(row.names(per_community))] <- per_community

  block_start <- 1

  while(block_start <= n_communities){

    block_end <- block_start
    block_work <- work[block_start]

    while(block_end < n_communities &&
          block_end - block_start + 1 < max_communities &&
          block_work + work[block_end + 1] <= max_work){

      block_end <- block_end + 1
      block_work <- block_work + work[block_end]

    }

    n_block <- block_end - block_start + 1

    if(last_row[block_end] >= first_row[block_start]){

      rows <- first_row[block_start]:last_row[block_end]

      local <- community[rows] - block_start + 1

      #How many of each community's tips lie below each edge.
      counts <- tabulate((rep(local, structure$path_length[tip[rows]]) - 1) * n_edges +
                           unlist(structure$path_edges[tip[rows]]),
                         nbins = n_block * n_edges)

      block_richness <- rep(richness[block_start:block_end], each = n_edges)

      #An edge is on the connecting subtree when some, but not all, are below it.
      included <- which(counts > 0 & counts < block_richness)

      if(length(included) > 0){

        included_community <- (included - 1) %/% n_edges + 1
        included_edge <- (included - 1) %% n_edges + 1

        included_lengths <- lengths[included_edge, , drop = FALSE]

        pooled <- rowsum(included_lengths, included_community, reorder = FALSE)

        pd[block_start - 1 + as.integer(row.names(pooled)), ] <- pooled

        if(!is.null(rarefy_to)){

          below <- counts[included]
          total <- block_richness[included]

          #A community smaller than the target cannot be rarefied to it.
          enough <- total >= rarefy_to

          log_all <- lchoose(total[enough], rarefy_to)

          all_from <- function(x){

            out <- numeric(length(x))
            possible <- x >= rarefy_to
            out[possible] <- exp(lchoose(x[possible], rarefy_to) - log_all[possible])
            out

          }

          kept <- 1 - all_from(below[enough]) - all_from(total[enough] - below[enough])

          if(any(enough)){

            pooled <- rowsum(included_lengths[enough, , drop = FALSE] * kept,
                             included_community[enough],
                             reorder = FALSE)

            rarefied[block_start - 1 + as.integer(row.names(pooled)), ] <- pooled

          }

        }

      }

    }

    block_start <- block_end + 1

  }

  list(richness = richness, pd = pd, rarefied = rarefied)

}

#' Phylogenetic diversity and its standardized index for many communities
#'
#' Internal helper.  Ties the pieces above together for a set of trees. Trees
#' that share the first tree's topology are handled in one pass; any that do not
#' are handled on their own.
#'
#' The conventions follow `PhyloMeasures::pd.query()` under its uniform null
#' model, which this replaces: a community of fewer than two tips has PD 0 and
#' index 0, and where the null standard deviation is zero the index is the
#' undivided difference from the mean.
#'
#' @param trees A list of phylogenies with the same tips.
#' @param community Integer community index of each occurrence.
#' @param species Species name of each occurrence. Names absent from a tree are
#'   ignored, and repeats within a community are counted once.
#' @param n_communities Number of communities.
#' @param rarefy_to Optional richness to rarefy every community to before
#'   standardizing, which removes the dependence of the index on richness
#'   (Sandel 2018). See Details.
#' @details Sandel (2018) rarefies each community to a common richness and
#'   averages the index over many random rarefactions. For a fixed target
#'   richness the index is a linear function of PD, so that average is
#'
#'       (expected PD of a random subset of the community - mean) / sd
#'
#'   with the mean and sd those of the target richness, and the expected PD of a
#'   random subset given exactly by the null mean formula applied to the
#'   community's own subtree. No replication is needed, and none is done.
#' @return A list: `richness`, from the first tree; `pd` and `pdi`, communities
#'   by trees matrices; and, with `rarefy_to`, `rarefied_pd` and `rarefied_pdi`,
#'   the same shape, NA for communities of fewer species than the target.
#' @noRd
phylo_diversity <- function(trees, community, species, n_communities, rarefy_to = NULL){

  pd <- matrix(data = NA_real_, nrow = n_communities, ncol = length(trees))
  pdi <- pd
  rarefied_pd <- pd
  rarefied_pdi <- pd

  first <- trees[[1]]

  shared <- vapply(X = trees,
                   FUN = function(tree_i){

                     identical(tree_i$edge, first$edge) &&
                       identical(tree_i$tip.label, first$tip.label)

                   },
                   FUN.VALUE = logical(1))

  #Trees arranged differently from the first are done one at a time.
  for(other in which(!shared)){

    alone <- phylo_diversity(trees = trees[other],
                             community = community,
                             species = species,
                             n_communities = n_communities,
                             rarefy_to = rarefy_to)

    pd[, other] <- alone$pd
    pdi[, other] <- alone$pdi

    if(!is.null(rarefy_to)){

      rarefied_pd[, other] <- alone$rarefied_pd
      rarefied_pdi[, other] <- alone$rarefied_pdi

    }

  }

  structure <- pd_tree_structure(first)

  lengths <- vapply(X = trees[shared],
                    FUN = function(tree_i) as.numeric(tree_i$edge.length),
                    FUN.VALUE = numeric(structure$n_edges))

  lengths <- matrix(data = lengths, nrow = structure$n_edges)

  #One record per species per community, for species that are in the tree.
  tip <- match(species, first$tip.label)

  usable <- !is.na(tip) & !duplicated(data.frame(community, tip))

  observed <- pd_of_communities(structure = structure,
                                lengths = lengths,
                                community = community[usable],
                                tip = tip[usable],
                                n_communities = n_communities,
                                rarefy_to = rarefy_to)

  sums <- pd_branch_sums(structure = structure, lengths = lengths)

  #Subtract the null mean and divide by the null sd, tree by tree. Where the sd
  #is zero the difference is left undivided.
  standardize <- function(values, richness){

    moments <- pd_null_moments(structure = structure, sums = sums, richness = richness)

    centred <- sweep(values, 2, moments$mean)

    varies <- moments$sd > 0

    centred[, varies] <- sweep(centred[, varies, drop = FALSE], 2, moments$sd[varies], "/")

    centred

  }

  index <- observed$pd

  for(richness in unique(observed$richness)){

    members <- which(observed$richness == richness)

    index[members, ] <- standardize(observed$pd[members, , drop = FALSE], richness)

  }

  pd[, shared] <- observed$pd
  pdi[, shared] <- index

  if(is.null(rarefy_to)){

    return(list(richness = observed$richness, pd = pd, pdi = pdi))

  }

  #Every community is compared with the same null, that of the target richness.
  rarefied_pd[, shared] <- observed$rarefied
  rarefied_pdi[, shared] <- standardize(observed$rarefied, rarefy_to)

  list(richness = observed$richness,
       pd = pd,
       pdi = pdi,
       rarefied_pd = rarefied_pd,
       rarefied_pdi = rarefied_pdi)

}
