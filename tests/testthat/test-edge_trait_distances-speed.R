# Speed of the vectorized branch-length assignment against the edge-at-a-time
# loop it replaced (issue #1), on both complete and incomplete trait data.
# Helpers live in helper-branch_lengths.R.

# The observed speedup is two to three orders of magnitude and grows with tree
# size; this threshold is set far below that so the test flags a return to
# quadratic behaviour rather than ordinary timing noise.
minimum_speedup <- 10

test_that("phylopars() fills in missing tip values, so the branch-length step sees dense input either way", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 60, seed = 7)

  complete <- sim$traits
  incomplete <- add_missing_values(complete)

  expect_true(any(is.na(incomplete)))

  anc_complete <- Rphylopars::phylopars(trait_data = complete, tree = sim$tree)$anc_recon
  anc_incomplete <- Rphylopars::phylopars(trait_data = incomplete, tree = sim$tree)$anc_recon

  # Both reconstructions are the same shape and neither contains NAs, which is
  # why the cost of the branch-length step does not depend on whether the input
  # traits were complete.
  expect_equal(dim(anc_incomplete), dim(anc_complete))
  expect_false(any(is.na(anc_complete)))
  expect_false(any(is.na(anc_incomplete)))

})

test_that("vectorized assignment beats the edge-at-a-time loop on complete data", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 200, seed = 8)

  anc_recon <- Rphylopars::phylopars(trait_data = sim$traits, tree = sim$tree)$anc_recon

  loop_time <- time_per_call(reference_branch_lengths(sim$tree, anc_recon))
  vectorized_time <- time_per_call(edge_trait_distances(sim$tree, anc_recon))

  expect_equal(edge_trait_distances(sim$tree, anc_recon),
               reference_branch_lengths(sim$tree, anc_recon),
               tolerance = 1e-10)

  expect_gt(loop_time/vectorized_time, minimum_speedup)

})

test_that("vectorized assignment beats the edge-at-a-time loop on incomplete data", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 200, seed = 8)

  traits <- add_missing_values(sim$traits)

  anc_recon <- Rphylopars::phylopars(trait_data = traits, tree = sim$tree)$anc_recon

  loop_time <- time_per_call(reference_branch_lengths(sim$tree, anc_recon))
  vectorized_time <- time_per_call(edge_trait_distances(sim$tree, anc_recon))

  expect_equal(edge_trait_distances(sim$tree, anc_recon),
               reference_branch_lengths(sim$tree, anc_recon),
               tolerance = 1e-10)

  expect_gt(loop_time/vectorized_time, minimum_speedup)

})

test_that("the vectorized assignment does not scale quadratically with tree size", {

  skip_on_cran()

  small <- simulate_traits(n_tips = 250, seed = 9)
  large <- simulate_traits(n_tips = 1000, seed = 9)

  anc_small <- Rphylopars::phylopars(trait_data = small$traits, tree = small$tree)$anc_recon
  anc_large <- Rphylopars::phylopars(trait_data = large$traits, tree = large$tree)$anc_recon

  small_time <- time_per_call(edge_trait_distances(small$tree, anc_small))
  large_time <- time_per_call(edge_trait_distances(large$tree, anc_large))

  # A 4x larger tree costs 16x under the old loop.  Allow generous headroom over
  # the linear expectation of 4x while still catching a return to quadratic.
  expect_lt(large_time/small_time, 10)

})
