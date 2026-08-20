# Reference implementation: the original edge-at-a-time loop, kept here so the
# vectorized branch-length assignment in scale_branches_multidimensional() can be
# checked against it.  See issue #1.

reference_branch_lengths <- function(tree, anc_recon, rate = FALSE){

  row.names(anc_recon)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)

  output_branches <- numeric(length(tree$edge.length))

  for(i in 1:length(tree$edge.length)){

    node_1 <- tree$edge[i,][1]
    node_2 <- tree$edge[i,][2]

    value_1 <- anc_recon[which(row.names(anc_recon) == node_1),]
    value_2 <- anc_recon[which(row.names(anc_recon) == node_2),]

    bl <- stats::dist(rbind(value_1, value_2))[1]

    if(rate){

      bl <- bl/tree$edge.length[i]

    }

    output_branches[i] <- bl

  }

  output_branches

}

# Simulate a small tree with several correlated traits.

simulate_traits <- function(n_tips = 30, n_traits = 4, seed = 1){

  set.seed(seed)

  tree <- ape::rcoal(n = n_tips)
  tree$tip.label <- paste0("sp", 1:n_tips)

  # Correlated Brownian motion: a random covariance matrix shared across traits.
  loadings <- matrix(stats::rnorm(n_traits * n_traits), nrow = n_traits)
  sigma <- crossprod(loadings) + diag(n_traits)

  trait_values <- phytools::sim.corrs(tree = tree, vcv = sigma)

  traits <- data.frame(species = row.names(trait_values),
                       trait_values,
                       row.names = NULL,
                       stringsAsFactors = FALSE)

  names(traits)[-1] <- paste0("trait", 1:n_traits)

  list(tree = tree, traits = traits)

}

test_that("vectorized branch lengths match the edge-at-a-time reference (complete data)", {

  sim <- simulate_traits()

  anc_recon <- Rphylopars::phylopars(trait_data = sim$traits,
                                     tree = sim$tree)$anc_recon

  for(rate in c(FALSE, TRUE)){

    scaled <- scale_branches_multidimensional(tree = sim$tree,
                                              traits = sim$traits,
                                              rate = rate)

    expected <- reference_branch_lengths(tree = sim$tree,
                                         anc_recon = anc_recon,
                                         rate = rate)

    expect_equal(as.numeric(scaled$edge.length), expected, tolerance = 1e-10)

  }

})

test_that("vectorized branch lengths match the reference when traits are missing", {

  sim <- simulate_traits(seed = 2)

  traits <- sim$traits

  # Scatter NAs through the trait columns, leaving every species at least one
  # observed trait.
  set.seed(3)
  for(j in 2:ncol(traits)){

    traits[sample(x = nrow(traits), size = 5), j] <- NA

  }

  anc_recon <- Rphylopars::phylopars(trait_data = traits,
                                     tree = sim$tree)$anc_recon

  for(rate in c(FALSE, TRUE)){

    scaled <- scale_branches_multidimensional(tree = sim$tree,
                                              traits = traits,
                                              rate = rate)

    expected <- reference_branch_lengths(tree = sim$tree,
                                         anc_recon = anc_recon,
                                         rate = rate)

    expect_equal(as.numeric(scaled$edge.length), expected, tolerance = 1e-10)

  }

})

test_that("output is a valid phylo object and return_traits still works", {

  sim <- simulate_traits(seed = 4)

  scaled <- scale_branches_multidimensional(tree = sim$tree,
                                            traits = sim$traits)

  expect_s3_class(scaled, "phylo")
  expect_true(is.numeric(scaled$edge.length))
  expect_null(dim(scaled$edge.length))
  expect_length(scaled$edge.length, nrow(sim$tree$edge))
  expect_true(all(is.finite(scaled$edge.length)))

  with_traits <- scale_branches_multidimensional(tree = sim$tree,
                                                 traits = sim$traits,
                                                 return_traits = TRUE)

  expect_named(with_traits, c("tree", "traits"))
  expect_equal(with_traits$tree$edge.length, scaled$edge.length)
  expect_equal(nrow(with_traits$traits), length(sim$tree$tip.label))
  expect_equal(row.names(with_traits$traits), sim$tree$tip.label)

})
