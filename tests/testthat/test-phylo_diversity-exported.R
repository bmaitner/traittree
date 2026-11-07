# The exported interface to the phylogenetic diversity calculations:
# phylo_diversity() and phylo_diversity_moments(). The arithmetic underneath is
# tested in test-phylo_diversity.R and test-rarefied_pdi.R; these tests are
# about what users pass in and get back.

reference <- readRDS(test_path("fixtures", "pd_reference.rds"))

reference_trees <- c("ultrametric", "ragged", "tiny", "ladder")

# the tidy form of a communities by species matrix
as_occurrences <- function(communities){

  occupied <- which(communities == 1, arr.ind = TRUE)

  data.frame(species = colnames(communities)[occupied[, 2]],
             community = occupied[, 1],
             stringsAsFactors = FALSE)

}

test_that("a community matrix gives the values PhyloMeasures gave", {

  for(name in reference_trees){

    case <- reference[[name]]

    result <- phylo_diversity(tree = case$tree, communities = case$communities)

    expect_s3_class(result, "data.frame")
    expect_named(result, c("community", "richness", "pd", "pdi"))
    expect_equal(nrow(result), nrow(case$communities))

    expect_equal(result$pd, case$pd, tolerance = 1e-10, info = name)
    expect_lt(max(abs(result$pdi - case$pdi)), 1e-6)
    expect_equal(result$richness, unname(rowSums(case$communities)), info = name)

  }

})

test_that("occurrences in tidy form give the same answer as the matrix", {

  case <- reference$ragged

  # drop the empty community, which a table of occurrences cannot express
  occupied <- rowSums(case$communities) > 0

  from_matrix <- phylo_diversity(case$tree, case$communities[occupied, ])
  from_tidy <- phylo_diversity(case$tree, as_occurrences(case$communities[occupied, ]))

  # tidy communities come back in the order they first appear, so line them up
  from_tidy <- from_tidy[order(from_tidy$community), ]

  expect_equal(from_tidy$richness, from_matrix$richness)
  expect_equal(from_tidy$pd, from_matrix$pd)
  expect_equal(from_tidy$pdi, from_matrix$pdi)

  # factors are read as names, not as their integer codes
  as_factors <- as_occurrences(case$communities[occupied, ])
  as_factors$species <- factor(as_factors$species)

  with_factors <- phylo_diversity(case$tree, as_factors)

  expect_equal(with_factors$pd[order(with_factors$community)], from_matrix$pd)

})

test_that("communities keep their names, in the order given", {

  tree <- reference$ultrametric$tree

  occurrences <- data.frame(species = tree$tip.label[c(1, 2, 3, 4, 5, 6, 7)],
                            site = c("zeta", "zeta", "alpha", "alpha", "alpha", "mid", "mid"),
                            stringsAsFactors = FALSE)

  result <- phylo_diversity(tree, occurrences)

  expect_equal(result$community, c("zeta", "alpha", "mid"))
  expect_equal(result$richness, c(2, 3, 2))

  # row names of a matrix are used the same way
  wide <- matrix(0, nrow = 2, ncol = 4,
                 dimnames = list(c("north", "south"), tree$tip.label[1:4]))
  wide["north", 1:2] <- 1
  wide["south", 2:4] <- 1

  expect_equal(phylo_diversity(tree, wide)$community, c("north", "south"))

  # and without them communities are numbered
  rownames(wide) <- NULL

  expect_equal(phylo_diversity(tree, wide)$community, 1:2)

})

test_that("a data frame of species columns and a table are read as matrices", {

  tree <- reference$ultrametric$tree
  case <- reference$ultrametric

  expect_equal(phylo_diversity(tree, as.data.frame(case$communities))$pd,
               case$pd, tolerance = 1e-10)

  # table() of occurrences, as in the example
  occurrences <- as_occurrences(case$communities[rowSums(case$communities) > 0, ])
  tabulated <- table(occurrences$community, occurrences$species)

  from_table <- phylo_diversity(tree, tabulated)
  from_occurrences <- phylo_diversity(tree, occurrences)

  expect_equal(from_table$pd[match(from_occurrences$community, from_table$community)],
               from_occurrences$pd)

})

test_that("any non-zero entry counts as present, so abundances are accepted", {

  case <- reference$ragged

  abundances <- case$communities * 7.5
  logical_form <- case$communities == 1

  presence <- phylo_diversity(case$tree, case$communities)

  expect_equal(phylo_diversity(case$tree, abundances), presence)
  expect_equal(phylo_diversity(case$tree, logical_form), presence)

})

test_that("standardize = FALSE returns PD alone", {

  case <- reference$ragged

  result <- phylo_diversity(case$tree, case$communities, standardize = FALSE)

  expect_named(result, c("community", "richness", "pd"))
  expect_equal(result$pd, case$pd, tolerance = 1e-10)

  rarefied <- phylo_diversity(case$tree, case$communities, standardize = FALSE, rarefy_to = 3)

  expect_named(rarefied, c("community", "richness", "pd", "pd_rarefied"))

})

test_that("rarefy_to adds rarefied PD and index, NA below the target", {

  case <- reference$ultrametric

  result <- phylo_diversity(case$tree, case$communities, rarefy_to = 5)

  expect_named(result, c("community", "richness", "pd", "pdi", "pd_rarefied", "pdi_rarefied"))
  expect_equal(attr(result, "rarefied_to"), 5)

  below <- result$richness < 5

  expect_true(all(is.na(result$pd_rarefied[below])))
  expect_true(all(is.na(result$pdi_rarefied[below])))
  expect_false(anyNA(result$pd_rarefied[!below]))
  expect_false(anyNA(result$pdi_rarefied[!below]))

  # at the target the rarefied and plain values are the same
  at_target <- result$richness == 5

  expect_gt(sum(at_target), 0)
  expect_equal(result$pd_rarefied[at_target], result$pd[at_target], tolerance = 1e-12)
  expect_equal(result$pdi_rarefied[at_target], result$pdi[at_target], tolerance = 1e-10)

  # above it, a random subset spans less of the tree than the whole community
  expect_true(all(result$pd_rarefied[result$richness > 5] < result$pd[result$richness > 5]))

  # the plain columns are untouched by asking for rarefaction
  expect_equal(result[, 1:4], phylo_diversity(case$tree, case$communities))

})

test_that("rarefy_to = \"min\" is the smallest richness of two or more", {

  case <- reference$ultrametric

  result <- phylo_diversity(case$tree, case$communities, rarefy_to = "min")

  richness <- rowSums(case$communities)

  expect_equal(attr(result, "rarefied_to"), min(richness[richness >= 2]))
  expect_equal(!is.na(result$pdi_rarefied), unname(richness >= 2))

})

test_that("the rarefied PD is the mean of the rarefaction curve of the community", {

  # Rarefied PD is the expected PD of a random subset of the community, which
  # is the null mean on the community's own subtree.
  tree <- reference$ultrametric$tree

  set.seed(61)
  species <- sample(tree$tip.label, 14)

  result <- phylo_diversity(tree,
                            data.frame(species = species, community = "a"),
                            rarefy_to = 6)

  own_subtree <- ape::keep.tip(tree, species)

  expect_equal(result$pd_rarefied,
               phylo_diversity_moments(own_subtree, 6)$mean,
               tolerance = 1e-10)

})

test_that("several phylogenies return matrices that match one at a time", {

  case <- reference$ultrametric

  set.seed(62)

  trees <- lapply(1:3, function(i){

    tree_i <- case$tree
    tree_i$edge.length <- stats::runif(nrow(tree_i$edge), 0, 2)
    tree_i

  })

  names(trees) <- c("first", "second", "third")
  class(trees) <- "multiPhylo"

  result <- phylo_diversity(trees, case$communities, rarefy_to = 3)

  expect_type(result, "list")
  expect_named(result, c("community", "richness", "pd", "pdi",
                         "pd_rarefied", "pdi_rarefied", "rarefied_to"))

  expect_equal(dim(result$pd), c(nrow(case$communities), 3))
  expect_equal(colnames(result$pdi), c("first", "second", "third"))
  expect_equal(result$rarefied_to, 3)

  for(i in 1:3){

    alone <- phylo_diversity(trees[[i]], case$communities, rarefy_to = 3)

    expect_equal(unname(result$pd[, i]), alone$pd, tolerance = 1e-12)
    expect_equal(unname(result$pdi[, i]), alone$pdi, tolerance = 1e-9)
    expect_equal(unname(result$pdi_rarefied[, i]), alone$pdi_rarefied, tolerance = 1e-9)

  }

  # a plain list works like a multiPhylo, and unnamed trees are numbered
  unnamed <- unname(unclass(trees))

  expect_equal(colnames(phylo_diversity(unnamed, case$communities)$pd), c("1", "2", "3"))

})

test_that("a multiPhylo with its tip labels compressed into an attribute still works", {

  # ape can store one copy of the labels on the multiPhylo instead of in each
  # tree. Unclassing such an object leaves every tree without tips.
  case <- reference$ultrametric

  set.seed(64)

  trees <- lapply(1:3, function(i){

    tree_i <- case$tree
    tree_i$edge.length <- stats::runif(nrow(tree_i$edge), 0, 2)
    tree_i

  })

  class(trees) <- "multiPhylo"
  compressed <- ape::.compressTipLabel(trees)

  expect_null(unclass(compressed)[[1]]$tip.label)

  expect_equal(phylo_diversity(compressed, case$communities),
               phylo_diversity(trees, case$communities))

  expect_equal(phylo_diversity_moments(compressed, c(2, 10, 30)),
               phylo_diversity_moments(trees, c(2, 10, 30)))

})

test_that("species that are not in the phylogeny are ignored with a warning", {

  tree <- reference$ultrametric$tree

  occurrences <- data.frame(species = c(tree$tip.label[1:3], "Not_a_tip", "Nor_this"),
                            community = c(1, 1, 1, 1, 2),
                            stringsAsFactors = FALSE)

  expect_warning(result <- phylo_diversity(tree, occurrences),
                 "2 species in communities are not tips")

  # the community holding only an unknown species is kept, with nothing in it
  expect_equal(result$richness, c(3, 0))
  expect_equal(result$pd[2], 0)

  clean <- phylo_diversity(tree, occurrences[1:3, ])

  expect_equal(result$pd[1], clean$pd)

  expect_warning(phylo_diversity(tree, occurrences[1:4, ]), "1 species in communities is not a tip")

})

test_that("unusable input is rejected with a message that says why", {

  tree <- reference$ultrametric$tree
  communities <- reference$ultrametric$communities

  expect_error(phylo_diversity("not a tree", communities), "must be a phylo object")
  expect_error(phylo_diversity(list(tree, "nor this"), communities), "must be a phylo object")

  no_lengths <- tree
  no_lengths$edge.length <- NULL

  expect_error(phylo_diversity(no_lengths, communities), "needs a length for every branch")

  negative <- tree
  negative$edge.length[1] <- -1

  expect_error(phylo_diversity(negative, communities), "must not be negative")

  fewer_tips <- ape::drop.tip(tree, tree$tip.label[1:5])

  expect_error(phylo_diversity(list(tree, fewer_tips), communities), "must have the same tips")

  unnamed <- unname(communities)

  expect_error(phylo_diversity(tree, unnamed), "must be named by species")

  with_gap <- communities
  with_gap[1, 1] <- NA

  expect_error(phylo_diversity(tree, with_gap), "missing values")
  expect_error(phylo_diversity(tree, -communities), "negative values")

  expect_error(phylo_diversity(tree, matrix("a", 2, 2, dimnames = list(NULL, c("x", "y")))),
               "numeric or logical table")

  expect_error(phylo_diversity(tree, data.frame(species = c("a", NA), community = c(1, 1))),
               "missing values")

  expect_error(phylo_diversity(tree, communities, standardize = NA), "must be TRUE or FALSE")

  for(bad in list(1, 2.5, "max", c(2, 3))){

    expect_error(phylo_diversity(tree, communities, rarefy_to = bad), "rarefy_to must be NULL")

  }

  expect_error(phylo_diversity(tree, communities, rarefy_to = 61),
               "larger than the number of species")

  expect_warning(phylo_diversity(tree, communities[2:3, ], rarefy_to = 30),
                 "No community has 30 or more species")

})

test_that("the moments match PhyloMeasures and behave at the ends", {

  for(name in reference_trees){

    case <- reference[[name]]
    n_tips <- length(case$tree$tip.label)

    moments <- phylo_diversity_moments(case$tree, richness = 0:n_tips)

    expect_s3_class(moments, "data.frame")
    expect_named(moments, c("richness", "mean", "sd"))

    expect_equal(moments$mean, case$moments[, 1], tolerance = 1e-10, info = name)
    expect_equal(moments$sd[-(n_tips + 1)], case$moments[-(n_tips + 1), 2],
                 tolerance = 1e-9, info = name)

    # no subtree below two species; one possible sample of everything
    expect_equal(moments$mean[1:2], c(0, 0))
    expect_equal(moments$sd[c(1, 2, n_tips + 1)], c(0, 0, 0))

    # drawing every species spans every branch
    expect_equal(moments$mean[n_tips + 1], sum(case$tree$edge.length), tolerance = 1e-10)

    # the rarefaction curve never falls as species are added
    expect_true(all(diff(moments$mean) >= -1e-12), info = name)

  }

})

test_that("the moments standardize PD into the index", {

  case <- reference$ragged

  diversity <- phylo_diversity(case$tree, case$communities)

  null <- phylo_diversity_moments(case$tree, diversity$richness)

  varies <- null$sd > 0

  expect_equal(((diversity$pd - null$mean)/null$sd)[varies], diversity$pdi[varies],
               tolerance = 1e-10)

  # richness is returned as asked, repeats and order included
  expect_equal(null$richness, diversity$richness)

})

test_that("moments for several phylogenies come back as matrices", {

  tree <- reference$ladder$tree

  set.seed(63)

  other <- tree
  other$edge.length <- stats::runif(nrow(tree$edge), 0, 1)

  # one sharing the layout, one with its edges in another order
  reordered <- ape::reorder.phylo(other, order = "postorder")

  result <- phylo_diversity_moments(list(a = tree, b = other, c = reordered), richness = c(2, 5, 12))

  expect_named(result, c("richness", "mean", "sd"))
  expect_equal(dim(result$mean), c(3, 3))
  expect_equal(colnames(result$sd), c("a", "b", "c"))

  expect_equal(unname(result$mean[, "a"]), phylo_diversity_moments(tree, c(2, 5, 12))$mean)
  expect_equal(unname(result$sd[, "b"]), phylo_diversity_moments(other, c(2, 5, 12))$sd)

  # the same tree laid out differently has the same moments
  expect_equal(result$mean[, "c"], result$mean[, "b"], tolerance = 1e-10)
  expect_equal(result$sd[, "c"], result$sd[, "b"], tolerance = 1e-9)

})

test_that("moments reject richness values that make no sense", {

  tree <- reference$tiny$tree

  expect_error(phylo_diversity_moments(tree, 6), "whole numbers between 0 and the number of tips")
  expect_error(phylo_diversity_moments(tree, -1), "whole numbers between 0")
  expect_error(phylo_diversity_moments(tree, 2.5), "whole numbers between 0")
  expect_error(phylo_diversity_moments(tree, NA), "whole numbers between 0")
  expect_error(phylo_diversity_moments(tree, numeric(0)), "whole numbers between 0")
  expect_error(phylo_diversity_moments("tree", 2), "must be a phylo object")

})
