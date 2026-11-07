# Phylogenetic diversity (PD) and its standardized index (PDI), computed in the
# package rather than through PhyloMeasures.
#
# Three independent lines of evidence are used: values frozen from PhyloMeasures
# 2.1 (fixtures/pd_reference.rds, made by dev/pdi_research/make_pd_fixtures.R),
# a direct evaluation of the published formulas over every pair of edges, and
# simulation.

reference <- readRDS(test_path("fixtures", "pd_reference.rds"))

reference_trees <- c("ultrametric", "ragged", "tiny", "ladder")

# The long form the functions take, from a communities by species 0/1 matrix.
long_form <- function(communities){

  occupied <- which(communities == 1, arr.ind = TRUE)

  list(community = unname(occupied[, 1]),
       species = colnames(communities)[occupied[, 2]],
       n = nrow(communities))

}

# The published formulas evaluated the slow way, over every pair of edges, with
# none of the pooling the package code relies on (Nipperess and Matsen 2013).
moments_by_pairs <- function(tree, richness){

  n_tips <- length(tree$tip.label)
  n_edges <- nrow(tree$edge)
  w <- tree$edge.length

  below <- matrix(0, n_edges, n_tips)

  for(tip in seq_len(n_tips)){

    path <- ape::nodepath(tree, from = n_tips + 1, to = tip)
    below[match(path[-1], tree$edge[, 2]), tip] <- 1

  }

  size <- rowSums(below)
  shared <- below %*% t(below)

  h <- function(x) ifelse(x < richness, 0, exp(lchoose(pmax(x, 0), richness) - lchoose(n_tips, richness)))

  size_e <- matrix(size, n_edges, n_edges)
  size_f <- t(size_e)

  both <- 1 - h(n_tips - size_e) - h(n_tips - size_f) + h(n_tips - size_e - size_f)

  lower <- pmin(size_e, size_f)
  upper <- pmax(size_e, size_f)
  nested <- shared == lower

  both[nested] <- (1 - h(upper) - h(n_tips - lower) + h(upper - lower))[nested]

  mean_pd <- sum(w * (1 - h(size) - h(n_tips - size)))

  c(mean = mean_pd, sd = sqrt(max(sum(outer(w, w) * both) - mean_pd^2, 0)))

}

test_that("PD matches PhyloMeasures on every reference tree", {

  for(name in reference_trees){

    case <- reference[[name]]
    input <- long_form(case$communities)

    result <- pd_by_community(trees = list(case$tree),
                              community = input$community,
                              species = input$species,
                              n_communities = input$n)

    expect_equal(result$pd[, 1], case$pd, tolerance = 1e-10, info = name)
    expect_equal(result$richness, unname(rowSums(case$communities)), info = name)

  }

})

test_that("the null mean and standard deviation match PhyloMeasures at every richness", {

  for(name in reference_trees){

    case <- reference[[name]]
    n_tips <- length(case$tree$tip.label)

    structure <- pd_tree_structure(case$tree)
    sums <- pd_branch_sums(structure, matrix(case$tree$edge.length, ncol = 1))

    ours <- t(vapply(0:n_tips,
                     function(r) unlist(pd_null_moments(structure, sums, r)),
                     numeric(2)))

    expect_equal(ours[, "mean"], case$moments[, 1], tolerance = 1e-10, info = name)

    # When every tip is drawn there is one possible sample and the standard
    # deviation is exactly zero. PhyloMeasures returns rounding noise there (of
    # order 1e-8), so that richness is compared with the truth instead.
    full <- n_tips + 1

    expect_equal(ours[-full, "sd"], case$moments[-full, 2], tolerance = 1e-9, info = name)
    expect_identical(ours[[full, "sd"]], 0, info = name)
    expect_lt(case$moments[full, 2], 1e-6)

  }

})

test_that("PDI matches PhyloMeasures, including its conventions at the edges", {

  for(name in reference_trees){

    case <- reference[[name]]
    input <- long_form(case$communities)

    result <- pd_by_community(trees = list(case$tree),
                              community = input$community,
                              species = input$species,
                              n_communities = input$n)

    # absolute tolerance: the reference carries ~1e-8 of noise where richness
    # equals the size of the tree
    expect_lt(max(abs(result$pdi[, 1] - case$pdi)), 1e-6)

    richness <- rowSums(case$communities)

    # Fewer than two species: no subtree, so PD and PDI are both zero.
    expect_true(all(result$pd[richness < 2, 1] == 0), info = name)
    expect_true(all(result$pdi[richness < 2, 1] == 0), info = name)

    # Every species: PD is the whole tree, the null has no spread, and the index
    # is the undivided difference from the mean, zero up to rounding.
    everything <- richness == length(case$tree$tip.label)

    expect_true(all(abs(result$pdi[everything, 1]) < 1e-10), info = name)

  }

})

test_that("the moments agree with the published formulas summed over every pair of edges", {

  # An evaluation that shares no code and none of the pooling with the package.
  for(name in c("ragged", "ladder", "tiny")){

    tree <- reference[[name]]$tree
    n_tips <- length(tree$tip.label)

    structure <- pd_tree_structure(tree)
    sums <- pd_branch_sums(structure, matrix(tree$edge.length, ncol = 1))

    for(richness in unique(c(2, 3, ceiling(n_tips/2), n_tips - 1))){

      ours <- pd_null_moments(structure, sums, richness)
      slow <- moments_by_pairs(tree, richness)

      expect_equal(ours$mean, slow[["mean"]], tolerance = 1e-10, info = paste(name, richness))
      expect_equal(ours$sd, slow[["sd"]], tolerance = 1e-8, info = paste(name, richness))

    }

  }

})

test_that("the moments agree with simulation", {

  skip_on_cran()

  tree <- reference$ragged$tree
  n_tips <- length(tree$tip.label)
  richness <- 9

  set.seed(41)
  n_draws <- 4000

  draws <- pd_by_community(trees = list(tree),
                           community = rep(seq_len(n_draws), each = richness),
                           species = unlist(lapply(seq_len(n_draws),
                                                   function(i) sample(tree$tip.label, richness))),
                           n_communities = n_draws)

  structure <- pd_tree_structure(tree)
  sums <- pd_branch_sums(structure, matrix(tree$edge.length, ncol = 1))
  exact <- pd_null_moments(structure, sums, richness)

  # the mean to within five standard errors, the sd to within 5%
  expect_lt(abs(mean(draws$pd) - exact$mean), 5 * exact$sd/sqrt(n_draws))
  expect_lt(abs(stats::sd(draws$pd) - exact$sd)/exact$sd, 0.05)

  # so the index of randomly assembled communities is centred on zero with unit spread
  expect_lt(abs(mean(draws$pdi)), 0.08)
  expect_lt(abs(stats::sd(draws$pdi) - 1), 0.05)

})

test_that("PD is the spanning subtree, without the path to the root", {

  tree <- ape::read.tree(text = "((a:1,b:2):10,(c:3,d:4):20);")

  measure <- function(species){

    pd_by_community(list(tree), rep(1, length(species)), species, 1)$pd[1, 1]

  }

  # Sisters: just their own two branches, not the 10 above them.
  expect_equal(measure(c("a", "b")), 3)
  expect_equal(measure(c("c", "d")), 7)

  # Across the root: both paths up to it.
  expect_equal(measure(c("a", "c")), 1 + 10 + 20 + 3)
  expect_equal(measure(c("a", "b", "c", "d")), sum(tree$edge.length))

})

test_that("trees sharing a topology are measured together and give the same answers as alone", {

  tree <- reference$ultrametric$tree
  input <- long_form(reference$ultrametric$communities)

  set.seed(42)
  variants <- lapply(1:5, function(i){

    tree_i <- tree
    tree_i$edge.length <- stats::runif(nrow(tree$edge), 0, 3)
    tree_i

  })

  together <- pd_by_community(variants, input$community, input$species, input$n)

  for(i in seq_along(variants)){

    alone <- pd_by_community(variants[i], input$community, input$species, input$n)

    expect_equal(together$pd[, i], alone$pd[, 1], tolerance = 1e-12)
    expect_equal(together$pdi[, i], alone$pdi[, 1], tolerance = 1e-9)

  }

})

test_that("a tree arranged differently from the first is still handled correctly", {

  tree <- reference$ragged$tree
  input <- long_form(reference$ragged$communities)

  # The same tree with its edges in another order, and one with its tips
  # relabelled in another order: neither shares the first tree's layout.
  reordered <- ape::reorder.phylo(tree, order = "postorder")
  rotated <- ape::rotateConstr(tree, rev(tree$tip.label))

  expect_false(identical(reordered$edge, tree$edge))

  result <- pd_by_community(list(tree, reordered, rotated),
                            input$community, input$species, input$n)

  expect_equal(result$pd[, 2], result$pd[, 1], tolerance = 1e-12)
  expect_equal(result$pd[, 3], result$pd[, 1], tolerance = 1e-12)
  expect_equal(result$pdi[, 2], result$pdi[, 1], tolerance = 1e-9)
  expect_equal(result$pdi[, 3], result$pdi[, 1], tolerance = 1e-9)

})

test_that("repeated records and species missing from the tree do not change the answer", {

  tree <- reference$ultrametric$tree

  species <- tree$tip.label[c(3, 17, 40, 41)]

  clean <- pd_by_community(list(tree), rep(1, 4), species, 1)

  messy <- pd_by_community(list(tree),
                           rep(1, 8),
                           c(species, species[1:2], "not_in_tree", "nor_this"),
                           1)

  expect_equal(messy$richness, 4)
  expect_equal(messy$pd, clean$pd)
  expect_equal(messy$pdi, clean$pdi)

})

test_that("communities with no usable records come back as zero, in place", {

  tree <- reference$ultrametric$tree

  # community 2 of 4 has no records at all; community 3 only an unknown species
  result <- pd_by_community(list(tree),
                            community = c(1, 1, 3, 4, 4, 4),
                            species = c(tree$tip.label[1:2], "unknown", tree$tip.label[5:7]),
                            n_communities = 4)

  expect_equal(result$richness, c(2, 0, 0, 3))
  expect_equal(result$pd[2:3, 1], c(0, 0))
  expect_equal(result$pdi[2:3, 1], c(0, 0))
  expect_true(all(result$pd[c(1, 4), 1] > 0))

})

test_that("working in blocks gives the same answer as working in one piece", {

  tree <- reference$ultrametric$tree
  structure <- pd_tree_structure(tree)
  lengths <- cbind(tree$edge.length, rev(tree$edge.length))

  set.seed(43)
  n_communities <- 57
  sizes <- sample(0:12, n_communities, replace = TRUE)

  community <- rep(seq_len(n_communities), sizes)
  tip <- unlist(lapply(sizes, function(k) sample(length(tree$tip.label), k)))

  whole <- pd_of_communities(structure, lengths, community, tip, n_communities)

  # one community per block, a handful per block, and a limit on work that
  # single communities exceed
  for(limits in list(c(1, Inf), c(5, Inf), c(Inf, 1), c(7, 40))){

    pieces <- pd_of_communities(structure, lengths, community, tip, n_communities,
                                max_communities = limits[1], max_work = limits[2])

    expect_equal(pieces$pd, whole$pd, info = paste(limits, collapse = "/"))
    expect_equal(pieces$richness, whole$richness)

  }

  # and the order the records arrive in is irrelevant
  shuffled <- sample(length(community))

  expect_equal(pd_of_communities(structure, lengths, community[shuffled], tip[shuffled],
                                 n_communities)$pd,
               whole$pd)

})

test_that("the tree structure accounts for every edge and every ancestor", {

  for(name in reference_trees){

    tree <- reference[[name]]$tree
    structure <- pd_tree_structure(tree)

    # tips below each edge, against ape's own count
    expect_equal(structure$clade_size,
                 ape::node.depth(tree)[tree$edge[, 2]],
                 info = name)

    # each tip's path ends at the root and visits each edge once
    for(tip in c(1, length(tree$tip.label))){

      path <- ape::nodepath(tree, from = length(tree$tip.label) + 1, to = tip)

      expect_setequal(structure$path_edges[[tip]], match(path[-1], tree$edge[, 2]))

    }

    # one ancestor pair for every step between an edge and the root, counted
    # independently as each edge's depth below the root
    unit <- tree
    unit$edge.length <- rep(1, nrow(tree$edge))

    expect_equal(length(structure$descendants),
                 sum(ape::node.depth.edgelength(unit)[tree$edge[, 2]] - 1),
                 info = name)

    expect_true(all(structure$clade_size[structure$descendants] <=
                      structure$clade_size[structure$ancestors]), info = name)

  }

})

test_that("a star tree reduces PD to species richness", {

  # With every tip hanging from the root on a branch of length one, PD is the
  # number of species, whatever they are, so the null has no spread at all.
  star <- ape::stree(12, type = "star")
  star$edge.length <- rep(1, 12)

  result <- pd_by_community(list(star),
                            community = rep(1:3, c(2, 5, 12)),
                            species = c(star$tip.label[1:2], star$tip.label[3:7], star$tip.label),
                            n_communities = 3)

  expect_equal(result$pd[, 1], c(2, 5, 12))

  # The variance is zero in truth and rounding noise in practice, so the index
  # is noise over noise: tiny, but not bit-exact zero.
  expect_true(all(abs(result$pdi[, 1]) < 1e-6))

})

test_that("a genuinely tiny variance is kept rather than rounded away", {

  # Four tips on a star, one branch longer by 2e-7: the null sd at richness 2
  # is 1e-7, a few parts in 1e7 of the mean, which is at the edge of what the
  # subtraction in the variance can resolve. It must survive, so that a pair
  # of the short branches scores about -1 and not about zero.
  star <- ape::stree(4, type = "star")
  star$tip.label <- letters[1:4]
  star$edge.length <- c(1, 1, 1, 1 + 2e-7)

  null <- phylo_diversity_moments(star, 2)

  expect_equal(null$mean, 2 + 1e-7, tolerance = 1e-12)
  expect_gt(null$sd, 0.5e-7)
  expect_lt(null$sd, 2e-7)

  pair <- phylo_diversity(star, data.frame(species = c("a", "b"), community = 1))

  expect_lt(pair$pdi, -0.5)
  expect_gt(pair$pdi, -2)

  # A negative value, which is pure rounding, never reaches the square root.
  equal <- star
  equal$edge.length <- rep(1, 4)

  expect_false(anyNA(phylo_diversity_moments(equal, 0:4)$sd))

})
