# percent = TRUE divides by the ancestral value at the start of each branch.
# Helpers live in helper-branch_lengths.R.

test_that("strictly positive, well separated denominators do not warn", {

  expect_no_warning(warn_unstable_percent_change(seq(from = 10, to = 20, length.out = 50)))

})

test_that("negative denominators warn", {

  values <- seq(from = -5, to = 5, length.out = 51)

  expect_warning(warn_unstable_percent_change(values), "are negative")
  expect_warning(warn_unstable_percent_change(values), "25 of 51")

})

test_that("denominators near zero warn even when all are positive", {

  # Spread of roughly 1, with three values pressed up against zero.
  values <- c(stats::runif(50, min = 1, max = 3), 1e-8, 1e-9, 1e-10)

  expect_warning(warn_unstable_percent_change(values), "of the trait's spread of zero")
  expect_warning(warn_unstable_percent_change(values), "3 of 53")

})

test_that("the near zero threshold is respected", {

  values <- c(rep(1, 20), rep(3, 20), 0.005)

  # the spread here is about 1.05, so 0.005 is inside 1% of it but not 0.1%.
  expect_warning(warn_unstable_percent_change(values, near_zero_fraction = 0.01), "1 of 41")
  expect_no_warning(warn_unstable_percent_change(values, near_zero_fraction = 0.001))

})

test_that("matrix input is assessed trait by trait", {

  # One well behaved trait, one straddling zero.
  denominator <- cbind(seq(from = 10, to = 20, length.out = 40),
                       seq(from = -4, to = 4, length.out = 40))

  expect_warning(warn_unstable_percent_change(denominator), "are negative")

  # Counts are over the whole matrix, so the well behaved column is included in
  # the total but contributes nothing to the count.
  expect_warning(warn_unstable_percent_change(denominator), "of 80")

})

test_that("degenerate input does not warn or error", {

  expect_no_warning(warn_unstable_percent_change(rep(5, 10)))
  expect_no_warning(warn_unstable_percent_change(1))

})

test_that("both scale_branches_by_traits functions warn when percent = TRUE", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, seed = 131)

  # Simulated Brownian traits are centred on zero, so ancestral values on both
  # sides of zero are the norm rather than an edge case.
  trait_vector <- sim$traits$trait1
  names(trait_vector) <- sim$traits$species

  expect_warning(scale_branches_by_traits_fastAnc(sim$tree, trait_vector, percent = TRUE),
                 "percent = TRUE divides by the ancestral value")

  expect_warning(suppressMessages(
    scale_branches_by_traits_rphylopars(sim$tree, sim$traits, percent = TRUE)),
    "percent = TRUE divides by the ancestral value")

})

test_that("percent = FALSE stays silent", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, seed = 131)

  trait_vector <- sim$traits$trait1
  names(trait_vector) <- sim$traits$species

  expect_no_warning(scale_branches_by_traits_fastAnc(sim$tree, trait_vector, percent = FALSE))

  expect_no_warning(suppressMessages(
    scale_branches_by_traits_rphylopars(sim$tree, sim$traits, percent = FALSE)))

})

test_that("strictly positive traits pass through percent = TRUE without warning", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, seed = 132)

  # Shift well clear of zero, as a strictly positive trait would be.
  positive <- sim$traits
  positive[, -1] <- positive[, -1] + 100

  trait_vector <- positive$trait1
  names(trait_vector) <- positive$species

  expect_no_warning(scale_branches_by_traits_fastAnc(sim$tree, trait_vector, percent = TRUE))

})
