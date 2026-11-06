# PhyloMeasures is a suggested dependency, since it was archived from CRAN and
# is only needed by the phylogenetic diversity functions.

test_that("the check passes when the package is available", {

  expect_true(require_phylomeasures(installed = TRUE))
  expect_invisible(require_phylomeasures(installed = TRUE))

})

test_that("the failure names the package and how to install it", {

  # `installed` is an argument so this path can be exercised without hiding the
  # package from the session.
  expect_error(require_phylomeasures(installed = FALSE), "PhyloMeasures is needed")
  expect_error(require_phylomeasures(installed = FALSE), "archived from CRAN")
  expect_error(require_phylomeasures(installed = FALSE), "install_url")
  expect_error(require_phylomeasures(installed = FALSE), "Archive/PhyloMeasures")

})

test_that("the real check agrees with requireNamespace", {

  available <- requireNamespace("PhyloMeasures", quietly = TRUE)

  if(available){

    expect_true(require_phylomeasures())

  }else{

    expect_error(require_phylomeasures(), "PhyloMeasures is needed")

  }

})

test_that("the rest of the package does not need PhyloMeasures", {

  # Nothing outside calculate_pd_metric_rasters() should touch it, so the trait
  # scaling path must not go anywhere near the guard.
  sim <- simulate_traits(n_tips = 20, seed = 201)

  expect_s3_class(scale_branches_multidimensional(sim$tree, sim$traits), "phylo")

  expect_s3_class(data_matching(phylogeny = sim$tree,
                                occurrences = data.frame(
                                  species = sim$tree$tip.label[1:10],
                                  cell = 1:10,
                                  stringsAsFactors = FALSE))$phylogeny,
                  "phylo")

})

test_that("bad arguments are reported before the dependency is checked", {

  # Argument errors are about the call, so they should not depend on whether an
  # optional package happens to be installed.
  sim <- simulate_traits(n_tips = 15, seed = 202)

  template <- raster::raster(nrows = 4, ncols = 4)

  occurrences <- data.frame(species = sim$tree$tip.label[1:10],
                            cell = rep(1:5, 2),
                            stringsAsFactors = FALSE)

  expect_error(calculate_pd_metric_rasters(occurrences, sim$tree, sim$traits,
                                           template, n_trees = 0),
               "n_trees must be a single positive number")

  expect_error(calculate_pd_metric_rasters(occurrences, sim$tree, sim$traits,
                                           template, probs = 0.5),
               "probs must be two probabilities")

  expect_error(calculate_pd_metric_rasters(occurrences, sim$tree, traits = NULL,
                                           template.raster = template),
               "Supply either traits")

})
