# Choosing whether to estimate a within-species (phenotypic) error term.
# Helpers live in helper-branch_lengths.R.

test_that("replication is detected per trait, not per row", {

  sim <- simulate_traits(n_tips = 20, seed = 5)

  # One row per species: nothing to estimate a within-species term from.
  expect_false(has_replicate_observations(sim$traits))

  # Several rows per species, but staggered NAs leave one observation of each
  # trait.  A rule counting rows per species would wrongly call this replicated.
  staggered <- sim$traits[rep(1:5, each = 2), ]
  staggered[seq(1, 10, by = 2), 2:3] <- NA
  staggered[seq(2, 10, by = 2), 4:5] <- NA

  expect_equal(max(table(staggered$species)), 2)
  expect_false(has_replicate_observations(staggered))

  # Genuine replication for a single species and a single trait is enough.
  one_replicate <- sim$traits
  extra_row <- one_replicate[1, ]
  extra_row[, 3:ncol(extra_row)] <- NA
  one_replicate <- rbind(one_replicate, extra_row)

  expect_true(has_replicate_observations(one_replicate))

})

test_that("has_replicate_observations handles degenerate input", {

  expect_false(has_replicate_observations(data.frame(species = character(0), trait1 = numeric(0))))
  expect_false(has_replicate_observations(data.frame(species = c("a", "a"))))
  expect_false(has_replicate_observations(data.frame(species = c("a", "a"), trait1 = c(NA, NA))))

})

test_that("resolve_pheno_error follows the data when pheno_error is NULL", {

  sim <- simulate_traits(n_tips = 20, seed = 5)

  replicated <- rbind(sim$traits, sim$traits[1:3, ])

  expect_false(resolve_pheno_error(traits = sim$traits, pheno_error = NULL))
  expect_true(resolve_pheno_error(traits = replicated, pheno_error = NULL))

  # NULL is the default.
  expect_false(resolve_pheno_error(traits = sim$traits))

})

test_that("an explicit pheno_error overrides what the data suggest", {

  sim <- simulate_traits(n_tips = 20, seed = 5)

  replicated <- rbind(sim$traits, sim$traits[1:3, ])

  expect_true(resolve_pheno_error(traits = sim$traits, pheno_error = TRUE))
  expect_false(resolve_pheno_error(traits = replicated, pheno_error = FALSE))

})

test_that("resolve_pheno_error rejects values that are not NULL, TRUE or FALSE", {

  sim <- simulate_traits(n_tips = 20, seed = 5)

  expect_error(resolve_pheno_error(sim$traits, pheno_error = NA), "must be NULL, TRUE, or FALSE")
  expect_error(resolve_pheno_error(sim$traits, pheno_error = "TRUE"), "must be NULL, TRUE, or FALSE")
  expect_error(resolve_pheno_error(sim$traits, pheno_error = c(TRUE, FALSE)), "must be NULL, TRUE, or FALSE")
  expect_error(resolve_pheno_error(sim$traits, pheno_error = 1), "must be NULL, TRUE, or FALSE")

})

# The value that actually reaches phylopars() is what matters, so intercept the
# call rather than inferring the setting from the results.

capture_pheno_error <- function(tree, traits, ...){

  captured <- NULL

  testthat::local_mocked_bindings(
    phylopars = function(trait_data, tree, pheno_error, ...){

      captured <<- pheno_error

      n_tips <- length(tree$tip.label)
      anc_recon <- matrix(0, nrow = 2 * n_tips - 1, ncol = ncol(trait_data) - 1)
      row.names(anc_recon) <- c(tree$tip.label, as.character((n_tips + 1):(2 * n_tips - 1)))

      list(anc_recon = anc_recon)

    },
    .package = "traittree"
  )

  scale_branches_multidimensional(tree = tree, traits = traits, ...)

  captured

}

test_that("scale_branches_multidimensional passes the resolved setting to phylopars", {

  sim <- simulate_traits(n_tips = 20, seed = 6)

  replicated <- rbind(sim$traits, sim$traits[1:3, ])

  # Chosen from the data by default...
  expect_false(capture_pheno_error(sim$tree, sim$traits))
  expect_true(capture_pheno_error(sim$tree, replicated))

  # ...and overridden when given explicitly.
  expect_true(capture_pheno_error(sim$tree, sim$traits, pheno_error = TRUE))
  expect_false(capture_pheno_error(sim$tree, replicated, pheno_error = FALSE))

})

test_that("species dropped from the tree do not drive the choice", {

  sim <- simulate_traits(n_tips = 20, seed = 6)

  # The only replication belongs to a species that is not in the tree, so it is
  # filtered out before the setting is chosen.
  absent <- sim$traits[1:2, ]
  absent$species <- "not_in_tree"

  expect_false(capture_pheno_error(sim$tree, rbind(sim$traits, absent)))

})

test_that("scale_branches_multidimensional runs end to end on replicated data", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 40, seed = 5)

  set.seed(12)
  extra <- do.call(rbind, replicate(3, sim$traits[1:10, ], simplify = FALSE))
  extra[, -1] <- extra[, -1] +
    matrix(stats::rnorm(nrow(extra) * (ncol(extra) - 1), sd = 0.1), nrow = nrow(extra))

  replicated <- rbind(sim$traits, extra)

  expect_true(has_replicate_observations(replicated))

  scaled <- scale_branches_multidimensional(tree = sim$tree, traits = replicated)

  expect_s3_class(scaled, "phylo")
  expect_length(scaled$edge.length, nrow(sim$tree$edge))
  expect_true(all(is.finite(scaled$edge.length)))

  with_traits <- scale_branches_multidimensional(tree = sim$tree,
                                                 traits = replicated,
                                                 return_traits = TRUE)

  expect_equal(nrow(with_traits$traits), length(sim$tree$tip.label))
  expect_equal(row.names(with_traits$traits), sim$tree$tip.label)

})

test_that("with one observation per species the new default reproduces the old results", {

  skip_on_cran()

  # phylopars() left pheno_error missing forces it to TRUE via its
  # pheno_correlated default, which is what this package used to do.  With a
  # single observation per species the term is not identifiable, so switching it
  # off must not move the results.
  for(seed in c(1, 2)){

    sim <- simulate_traits(seed = seed)

    for(traits in list(sim$traits, add_missing_values(sim$traits, proportion = 1/6))){

      old_default <- Rphylopars::phylopars(trait_data = traits, tree = sim$tree)$anc_recon
      new_default <- Rphylopars::phylopars(trait_data = traits, tree = sim$tree,
                                           pheno_error = FALSE)$anc_recon

      expect_equal(edge_trait_distances(sim$tree, new_default),
                   edge_trait_distances(sim$tree, old_default),
                   tolerance = 1e-10)

    }

  }

})

test_that("sparse replication warns, since one duplicated row would otherwise flip the model silently", {

  sim <- simulate_traits(n_tips = 100, seed = 11)

  # A single species out of 100 with a repeat measurement: enough to switch the
  # model specification, sparse enough to look like a data-entry artifact.
  barely <- rbind(sim$traits, sim$traits[1, ])

  expect_warning(result <- resolve_pheno_error(traits = barely),
                 "sparse replication")
  expect_true(result)

  # The warning names the counts, so the user can judge whether it was intended.
  expect_warning(resolve_pheno_error(traits = barely), "1 of 100 species")

})

test_that("ample replication does not warn", {

  sim <- simulate_traits(n_tips = 100, seed = 11)

  ample <- rbind(sim$traits, sim$traits[1:40, ])

  expect_no_warning(result <- resolve_pheno_error(traits = ample))
  expect_true(result)

})

test_that("no warning when the user sets pheno_error explicitly", {

  sim <- simulate_traits(n_tips = 100, seed = 11)

  barely <- rbind(sim$traits, sim$traits[1, ])

  # Having chosen the setting, the user does not need to be told about it.
  expect_no_warning(expect_true(resolve_pheno_error(barely, pheno_error = TRUE)))
  expect_no_warning(expect_false(resolve_pheno_error(barely, pheno_error = FALSE)))

})

test_that("the sparse replication threshold is respected at its boundary", {

  sim <- simulate_traits(n_tips = 100, seed = 11)

  # Exactly 5 of 100 species replicated: at the default threshold, not below it.
  at_threshold <- rbind(sim$traits, sim$traits[1:5, ])

  expect_no_warning(resolve_pheno_error(at_threshold))

  expect_warning(resolve_pheno_error(at_threshold, sparse_replication_threshold = 0.1),
                 "5 of 100 species")

  expect_no_warning(resolve_pheno_error(at_threshold, sparse_replication_threshold = 0))

})

test_that("the warning reaches callers of scale_branches_multidimensional", {

  sim <- simulate_traits(n_tips = 100, seed = 11)

  barely <- rbind(sim$traits, sim$traits[1, ])

  expect_warning(captured <- capture_pheno_error(sim$tree, barely), "sparse replication")
  expect_true(captured)

})
