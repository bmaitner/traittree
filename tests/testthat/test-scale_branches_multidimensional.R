# Correctness of the vectorized branch-length assignment (issue #1).
# Helpers live in helper-branch_lengths.R.

test_that("vectorized branch lengths match the edge-at-a-time reference (complete data)", {

  sim <- simulate_traits()

  # Match the setting the package resolves for this data, so the comparison is
  # about the branch-length assignment and nothing else.
  anc_recon <- Rphylopars::phylopars(trait_data = sim$traits,
                                     tree = sim$tree,
                                     pheno_error = resolve_pheno_error(sim$traits))$anc_recon

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

  traits <- add_missing_values(sim$traits, proportion = 1/6)

  anc_recon <- Rphylopars::phylopars(trait_data = traits,
                                     tree = sim$tree,
                                     pheno_error = resolve_pheno_error(traits))$anc_recon

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
