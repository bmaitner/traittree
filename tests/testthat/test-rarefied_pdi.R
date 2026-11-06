# The rarefied phylogenetic diversity index: PDI with its dependence on species
# richness removed, following Sandel (2018), but with the average over
# rarefactions computed exactly instead of by replication.

reference <- readRDS(test_path("fixtures", "pd_reference.rds"))

test_that("the closed form equals the average over every possible rarefaction", {

  # Sandel averages the index over many random rarefactions. Enumerating all of
  # them, rather than sampling, gives the quantity that estimates, exactly.
  tree <- reference$ultrametric$tree

  set.seed(51)

  for(case in list(c(5, 3), c(8, 3), c(9, 4), c(7, 2))){

    richness <- case[1]
    target <- case[2]

    species <- sample(tree$tip.label, richness)

    ours <- phylo_diversity(list(tree), rep(1, richness), species, 1, rarefy_to = target)

    subsets <- utils::combn(species, target)

    every_rarefaction <- phylo_diversity(list(tree),
                                         rep(seq_len(ncol(subsets)), each = target),
                                         as.vector(subsets),
                                         ncol(subsets))

    expect_equal(ours$rarefied_pdi[1, 1], mean(every_rarefaction$pdi[, 1]),
                 tolerance = 1e-10, info = paste(case, collapse = " to "))

    expect_equal(ours$rarefied_pd[1, 1], mean(every_rarefaction$pd[, 1]),
                 tolerance = 1e-10, info = paste(case, collapse = " to "))

  }

})

test_that("the same holds on a tree with polytomies and zero-length branches", {

  tree <- reference$ragged$tree

  set.seed(52)
  species <- sample(tree$tip.label, 8)

  ours <- phylo_diversity(list(tree), rep(1, 8), species, 1, rarefy_to = 4)

  subsets <- utils::combn(species, 4)

  every_rarefaction <- phylo_diversity(list(tree),
                                       rep(seq_len(ncol(subsets)), each = 4),
                                       as.vector(subsets),
                                       ncol(subsets))

  expect_equal(ours$rarefied_pdi[1, 1], mean(every_rarefaction$pdi[, 1]), tolerance = 1e-10)

})

test_that("a community already at the target is unchanged, and a smaller one is NA", {

  tree <- reference$ultrametric$tree

  result <- phylo_diversity(list(tree),
                            community = rep(1:3, c(6, 3, 12)),
                            species = c(tree$tip.label[1:6], tree$tip.label[10:12], tree$tip.label[20:31]),
                            n_communities = 3,
                            rarefy_to = 6)

  # at the target: there is one possible rarefaction, the community itself
  expect_equal(result$rarefied_pdi[1, 1], result$pdi[1, 1], tolerance = 1e-12)
  expect_equal(result$rarefied_pd[1, 1], result$pd[1, 1], tolerance = 1e-12)

  # below it: nothing to rarefy to
  expect_true(is.na(result$rarefied_pdi[2, 1]))
  expect_true(is.na(result$rarefied_pd[2, 1]))

  # above it: a real value, and a smaller expected PD than the whole community
  expect_false(is.na(result$rarefied_pdi[3, 1]))
  expect_lt(result$rarefied_pd[3, 1], result$pd[3, 1])

  # the plain values do not depend on whether rarefaction was requested
  plain <- phylo_diversity(list(tree),
                           community = rep(1:3, c(6, 3, 12)),
                           species = c(tree$tip.label[1:6], tree$tip.label[10:12], tree$tip.label[20:31]),
                           n_communities = 3)

  expect_equal(result$pd, plain$pd)
  expect_equal(result$pdi, plain$pdi)

})

test_that("rarefied values agree whether trees are measured together or apart", {

  tree <- reference$ultrametric$tree

  set.seed(53)

  variants <- lapply(1:4, function(i){

    tree_i <- tree
    tree_i$edge.length <- stats::runif(nrow(tree$edge), 0, 2)
    tree_i

  })

  sizes <- c(4, 9, 15, 30)
  community <- rep(seq_along(sizes), sizes)
  species <- unlist(lapply(sizes, function(k) sample(tree$tip.label, k)))

  together <- phylo_diversity(variants, community, species, length(sizes), rarefy_to = 5)

  # a tree laid out differently goes down the one-at-a-time path
  reordered <- ape::reorder.phylo(variants[[2]], order = "postorder")

  mixed <- phylo_diversity(list(variants[[1]], reordered), community, species,
                           length(sizes), rarefy_to = 5)

  expect_equal(mixed$rarefied_pdi[, 1], together$rarefied_pdi[, 1], tolerance = 1e-9)
  expect_equal(mixed$rarefied_pdi[, 2], together$rarefied_pdi[, 2], tolerance = 1e-9)

  expect_true(is.na(together$rarefied_pdi[1, 1]))
  expect_false(anyNA(together$rarefied_pdi[-1, ]))

})

test_that("rarefaction removes the richness-dependence Sandel (2018) describes", {

  skip_on_cran()

  # One filter, favouring a clade of roughly 40% of the tips fivefold, applied
  # to communities that differ only in how many species they hold. Richness is
  # kept well below the size of the favoured clade, as in the paper: beyond it
  # the filter saturates and the comparison means something else.
  #
  # On any one tree the plain index is erratic, so the comparison is pooled over
  # several: how far the mean index moves across richness levels, plain against
  # rarefied. Across seeds the rarefied index moves 20 to 30% as far.
  set.seed(2018)

  n_tips <- 400
  levels <- c(10, 25, 50)

  movement <- vapply(X = 1:6,
                     FUN = function(i){

                       tree <- ape::rcoal(n_tips)
                       tree$tip.label <- paste0("sp", seq_len(n_tips))

                       clade_sizes <- ape::node.depth(tree)
                       internal <- (n_tips + 1):(2 * n_tips - 1)
                       node <- internal[which.min(abs(clade_sizes[internal] - 0.4 * n_tips))]

                       favoured <- ape::extract.clade(tree, node)$tip.label
                       suitability <- ifelse(tree$tip.label %in% favoured, 1, 0.2)

                       richness <- rep(levels, each = 40)

                       community <- rep(seq_along(richness), richness)
                       species <- unlist(lapply(richness,
                                                function(k) sample(tree$tip.label, k, prob = suitability)))

                       result <- phylo_diversity(list(tree), community, species,
                                                 length(richness), rarefy_to = min(levels))

                       plain <- tapply(result$pdi[, 1], richness, mean)
                       rarefied <- tapply(result$rarefied_pdi[, 1], richness, mean)

                       c(plain = diff(range(plain)), rarefied = diff(range(rarefied)))

                     },
                     FUN.VALUE = numeric(2))

  movement <- rowMeans(movement)

  # The plain index shifts substantially with richness under the same filter...
  expect_gt(movement[["plain"]], 0.5)

  # ...and the rarefied index by well under half as much.
  expect_lt(movement[["rarefied"]], 0.5 * movement[["plain"]])

})

test_that("rarefy_to adds rarefied layers to the rasters and leaves the others alone", {

  template <- raster::raster(nrows = 10, ncols = 10)

  arguments <- list(occurrences = example_occurrences,
                    phylogeny = example_tree,
                    traits = NULL,
                    template.raster = template,
                    trait_phylo = reference$pipeline$trait_tree)

  plain <- do.call(calculate_pd_metric_rasters, arguments)
  rarefied <- do.call(calculate_pd_metric_rasters, c(arguments, list(rarefy_to = 4)))

  expect_equal(names(rarefied$pd_stack),
               c("pd_time", "pdi_time", "pd_traits", "pdi_traits",
                 "pdi_time_rarefied", "pdi_traits_rarefied"))

  expect_equal(rarefied$rarefied_to, 4)
  expect_null(plain$rarefied_to)

  values <- raster::values(rarefied$pd_stack)

  expect_equal(values[, 1:4], raster::values(plain$pd_stack))

  # Exactly the cells with at least four species carry a rarefied value.
  species_per_cell <- raster::values(richness_raster(template, example_occurrences))
  species_per_cell[is.na(species_per_cell)] <- 0

  expect_equal(!is.na(values[, "pdi_time_rarefied"]), species_per_cell >= 4)
  expect_equal(!is.na(values[, "pdi_traits_rarefied"]), species_per_cell >= 4)

  # In a cell of exactly four species the rarefied and plain indices coincide.
  at_target <- species_per_cell == 4

  expect_gt(sum(at_target), 0)
  expect_equal(values[at_target, "pdi_time_rarefied"], values[at_target, "pdi_time"],
               tolerance = 1e-10)

})

test_that("rarefied layers are summarised across a sample of trees like the others", {

  template <- raster::raster(nrows = 10, ncols = 10)

  result <- calculate_pd_metric_rasters(occurrences = example_occurrences,
                                        phylogeny = example_tree,
                                        traits = NULL,
                                        template.raster = template,
                                        trait_phylo = reference$pipeline$trait_sample,
                                        rarefy_to = 5,
                                        probs = c(0.1, 0.9))

  expect_true(all(c("pdi_time_rarefied",
                    "pdi_traits_rarefied_mean", "pdi_traits_rarefied_sd",
                    "pdi_traits_rarefied_q10", "pdi_traits_rarefied_q90") %in%
                    names(result$pd_stack)))

  values <- raster::values(result$pd_stack)
  present <- !is.na(values[, "pdi_traits_rarefied_mean"])

  expect_gt(sum(present), 0)
  expect_true(all(values[present, "pdi_traits_rarefied_q10"] <=
                    values[present, "pdi_traits_rarefied_q90"]))

  # the same cells are missing from every rarefied layer
  expect_equal(is.na(values[, "pdi_traits_rarefied_sd"]), !present)
  expect_equal(is.na(values[, "pdi_time_rarefied"]), !present)

})

test_that("rarefy_to = \"min\" uses the smallest richness that has a subtree", {

  template <- raster::raster(nrows = 10, ncols = 10)

  result <- calculate_pd_metric_rasters(occurrences = example_occurrences,
                                        phylogeny = example_tree,
                                        traits = NULL,
                                        template.raster = template,
                                        trait_phylo = reference$pipeline$trait_tree,
                                        rarefy_to = "min")

  species_per_cell <- raster::values(richness_raster(template, example_occurrences))

  expect_equal(result$rarefied_to, min(species_per_cell[species_per_cell >= 2], na.rm = TRUE))

  # every reported cell with a subtree gets a value
  values <- raster::values(result$pd_stack)

  expect_equal(!is.na(values[, "pdi_time_rarefied"]),
               !is.na(species_per_cell) & species_per_cell >= result$rarefied_to)

})

test_that("unusable values of rarefy_to are rejected or flagged", {

  template <- raster::raster(nrows = 10, ncols = 10)

  call_with <- function(rarefy_to){

    calculate_pd_metric_rasters(occurrences = example_occurrences,
                                phylogeny = example_tree,
                                traits = NULL,
                                template.raster = template,
                                trait_phylo = reference$pipeline$trait_tree,
                                rarefy_to = rarefy_to)

  }

  for(bad in list(1, 2.5, -3, NA_real_, c(3, 4), "max", TRUE)){

    expect_error(call_with(bad), "rarefy_to must be NULL")

  }

  # more species than the phylogeny holds
  expect_error(call_with(500), "larger than the number of species")

  # a legal target that no cell reaches
  expect_warning(result <- call_with(50), "No cell has 50 or more species")
  expect_true(all(is.na(raster::values(result$pd_stack[["pdi_time_rarefied"]]))))

})
