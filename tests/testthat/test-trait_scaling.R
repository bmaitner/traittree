# Traits enter the branch-length distance in whatever units they were measured
# in, unless they are scaled first.  Helpers live in helper-branch_lengths.R.

test_that("scaling makes branch lengths invariant to the units traits are measured in", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 40, n_traits = 3, seed = 141)

  # the same data with one trait in different units, e.g. kg rather than g
  rescaled <- sim$traits
  rescaled$trait1 <- rescaled$trait1 * 1000

  original <- scale_branches_multidimensional(sim$tree, sim$traits, scale_traits = TRUE)
  converted <- suppressWarnings(
    scale_branches_multidimensional(sim$tree, rescaled, scale_traits = TRUE))

  expect_equal(original$edge.length, converted$edge.length, tolerance = 1e-8)

  # Without scaling the unit change moves the answer, which is the reason the
  # default is TRUE.
  unscaled_original <- scale_branches_multidimensional(sim$tree, sim$traits, scale_traits = FALSE)
  unscaled_converted <- suppressWarnings(
    scale_branches_multidimensional(sim$tree, rescaled, scale_traits = FALSE))

  expect_lt(stats::cor(unscaled_original$edge.length, unscaled_converted$edge.length), 0.95)

})

test_that("scaling is the default", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, n_traits = 3, seed = 142)

  default <- scale_branches_multidimensional(sim$tree, sim$traits)
  explicit <- scale_branches_multidimensional(sim$tree, sim$traits, scale_traits = TRUE)

  expect_equal(default$edge.length, explicit$edge.length)

})

test_that("scaling the reconstruction matches rescaling the trait data", {

  skip_on_cran()

  # Brownian reconstruction is equivariant under a diagonal rescaling of the
  # traits, which is why the multipliers can be applied after the fit rather
  # than before it.
  sim <- simulate_traits(n_tips = 40, n_traits = 3, seed = 143)

  deviations <- trait_standard_deviations(sim$traits)

  pre_scaled <- sim$traits
  pre_scaled[, -1] <- sweep(as.matrix(sim$traits[, -1]), 2, deviations, "/")

  expect_equal(
    scale_branches_multidimensional(sim$tree, sim$traits, scale_traits = TRUE)$edge.length,
    scale_branches_multidimensional(sim$tree, pre_scaled, scale_traits = FALSE)$edge.length,
    tolerance = 1e-6)

})

test_that("return_traits stays in the units the caller supplied", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, n_traits = 3, seed = 144)

  result <- scale_branches_multidimensional(sim$tree, sim$traits,
                                            scale_traits = TRUE,
                                            return_traits = TRUE)

  observed <- as.matrix(sim$traits[match(sim$tree$tip.label, sim$traits$species), -1])

  expect_equal(unname(result$traits), unname(observed), tolerance = 1e-8)

})

test_that("weights are applied to the squared differences", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, n_traits = 3, seed = 145)

  equal <- scale_branches_multidimensional(sim$tree, sim$traits, scale_traits = TRUE)

  # A uniform weight is a constant factor on every branch, since the distance is
  # the square root of a weighted sum of squares.
  uniform <- scale_branches_multidimensional(sim$tree, sim$traits,
                                             scale_traits = TRUE,
                                             weights = rep(4, 3))

  expect_equal(uniform$edge.length, equal$edge.length * 2, tolerance = 1e-8)

  # Zero weight on a trait must give the same answer as dropping it.
  dropped <- scale_branches_multidimensional(sim$tree, sim$traits[, c("species", "trait1", "trait2")],
                                             scale_traits = TRUE)

  zeroed <- scale_branches_multidimensional(sim$tree, sim$traits,
                                            scale_traits = TRUE,
                                            weights = c(1, 1, 0))

  expect_equal(zeroed$edge.length, dropped$edge.length, tolerance = 1e-6)

})

test_that("weights are validated", {

  sim <- simulate_traits(n_tips = 15, n_traits = 3, seed = 146)

  expect_error(trait_scaling_multipliers(sim$traits, weights = c(1, 1)), "one per trait")
  expect_error(trait_scaling_multipliers(sim$traits, weights = c(1, 1, -1)), "non-negative")
  expect_error(trait_scaling_multipliers(sim$traits, weights = c(1, 1, NA)), "one per trait")
  expect_error(trait_scaling_multipliers(sim$traits, weights = c(1, 1, Inf)), "one per trait")
  expect_error(trait_scaling_multipliers(sim$traits, weights = c(0, 0, 0)), "all zero")
  expect_error(trait_scaling_multipliers(sim$traits, weights = c("1", "1", "1")), "numeric")

  expect_error(trait_scaling_multipliers(sim$traits, scale_traits = NA), "must be TRUE or FALSE")
  expect_error(trait_scaling_multipliers(sim$traits, scale_traits = "yes"), "must be TRUE or FALSE")

})

test_that("a trait with no variation is left alone rather than divided by zero", {

  sim <- simulate_traits(n_tips = 20, n_traits = 3, seed = 147)

  constant <- sim$traits
  constant$trait2 <- 7

  expect_warning(multipliers <- trait_scaling_multipliers(constant, scale_traits = TRUE),
                 "no variation across species")

  expect_true(all(is.finite(multipliers)))
  expect_equal(multipliers[2], 1)

})

test_that("standard deviations are taken over species, not observations", {

  sim <- simulate_traits(n_tips = 20, n_traits = 2, seed = 148)

  # Repeating a handful of species many times must not change the spread, since
  # replicates are averaged within a species first.
  replicated <- rbind(sim$traits,
                      do.call(rbind, replicate(20, sim$traits[1:3, ], simplify = FALSE)))

  expect_equal(trait_standard_deviations(replicated),
               trait_standard_deviations(sim$traits),
               tolerance = 1e-10)

})

test_that("using traits in their original units warns when their spreads differ wildly", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, n_traits = 3, seed = 149)

  lopsided <- sim$traits
  lopsided$trait1 <- lopsided$trait1 * 1000

  expect_warning(scale_branches_multidimensional(sim$tree, lopsided, scale_traits = FALSE),
                 "differ in spread by a factor of")

  # Comparable traits pass without comment...
  expect_no_warning(scale_branches_multidimensional(sim$tree, sim$traits, scale_traits = FALSE))

  # ...and so does anything where the caller has already made the choice.
  expect_no_warning(scale_branches_multidimensional(sim$tree, lopsided, scale_traits = TRUE))
  expect_no_warning(scale_branches_multidimensional(sim$tree, lopsided,
                                                    scale_traits = FALSE,
                                                    weights = c(1, 1, 1)))

})

test_that("with_variation scales its traits the same way", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, n_traits = 3, seed = 150)

  rescaled <- sim$traits
  rescaled$trait1 <- rescaled$trait1 * 1000

  set.seed(11)
  original <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, sim$traits, n_trees = 3))

  set.seed(11)
  converted <- suppressWarnings(suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, rescaled, n_trees = 3)))

  for(i in seq_along(original)){

    expect_equal(original[[i]]$edge.length, converted[[i]]$edge.length, tolerance = 1e-6)

  }

})
