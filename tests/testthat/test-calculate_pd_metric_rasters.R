# The batched rewrite of calculate_pd_metric_rasters(), checked against the
# original one-cell-at-a-time implementation in
# helper-reference-implementations.R.

# Build a small but non-trivial occurrence dataset over a template raster.
simulate_occurrences <- function(n_tips = 40, n_cells = 60, n_records = 300, seed = 111){

  sim <- simulate_traits(n_tips = n_tips, seed = seed)

  set.seed(seed)

  occurrences <- data.frame(species = sample(sim$tree$tip.label, n_records, replace = TRUE),
                            cell = sample(n_cells, n_records, replace = TRUE),
                            stringsAsFactors = FALSE)

  list(tree = sim$tree, traits = sim$traits, occurrences = occurrences)

}

test_that("batched pd/pdi rasters match the original cell-by-cell implementation", {

  skip_on_cran()

  fixture <- simulate_occurrences()

  template <- raster::raster(nrows = 10, ncols = 10)

  batched <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = fixture$traits,
                                template.raster = template))

  reference <- suppressMessages(
    reference_calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                          phylogeny = fixture$tree,
                                          traits = fixture$traits,
                                          template.raster = template))

  expect_equal(names(batched), names(reference))
  expect_equal(names(batched$pd_stack), names(reference$pd_stack))

  for(layer in names(batched$pd_stack)){

    expect_equal(raster::values(batched$pd_stack[[layer]]),
                 raster::values(reference$pd_stack[[layer]]),
                 tolerance = 1e-8,
                 info = layer)

  }

  # The phylogenies returned alongside the rasters should be untouched by the
  # change to how cells are queried.
  for(component in c("time_phylo", "trait_phylo", "rate_phylo")){

    expect_equal(batched[[component]]$tip.label, reference[[component]]$tip.label, info = component)
    expect_equal(batched[[component]]$edge.length, reference[[component]]$edge.length,
                 tolerance = 1e-8, info = component)

  }

})

test_that("cells with fewer than two records stay NA, and empty cells stay NA", {

  skip_on_cran()

  fixture <- simulate_occurrences(n_records = 80, n_cells = 60)

  template <- raster::raster(nrows = 10, ncols = 10)

  result <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = fixture$traits,
                                template.raster = template))

  values <- raster::values(result$pd_stack[["pd_time"]])

  # Recreate the occurrence filtering the function performs internally.
  matched <- data_matching(phylogeny = fixture$tree, occurrences = fixture$occurrences)
  records_per_cell <- table(as.character(matched$occurrences[, 2]))

  occupied <- as.numeric(names(records_per_cell))
  single_record <- occupied[records_per_cell < 2]
  multi_record <- occupied[records_per_cell >= 2]

  expect_true(all(is.na(values[single_record])))
  expect_true(all(!is.na(values[multi_record])))

  # cells with no occurrences at all
  expect_true(all(is.na(values[setdiff(seq_len(raster::ncell(template)), occupied)])))

})

test_that("verbose reports progress and is silent by default", {

  skip_on_cran()

  fixture <- simulate_occurrences(n_tips = 20, n_cells = 20, n_records = 100)

  template <- raster::raster(nrows = 5, ncols = 5)

  arguments <- list(occurrences = fixture$occurrences,
                    phylogeny = fixture$tree,
                    traits = fixture$traits,
                    template.raster = template)

  # scale_branches_by_traits_rphylopars is not involved here, but the
  # reconstruction still emits nothing, so the only messages are ours.
  expect_message(do.call(calculate_pd_metric_rasters, c(arguments, list(verbose = TRUE))),
                 "cells")

  expect_no_message(do.call(calculate_pd_metric_rasters, arguments))

})

test_that("n_trees > 1 replaces the trait layers with a summary across trees", {

  skip_on_cran()

  fixture <- simulate_occurrences()

  template <- raster::raster(nrows = 10, ncols = 10)

  arguments <- list(occurrences = fixture$occurrences,
                    phylogeny = fixture$tree,
                    traits = fixture$traits,
                    template.raster = template)

  point <- suppressMessages(do.call(calculate_pd_metric_rasters, arguments))
  sampled <- suppressMessages(do.call(calculate_pd_metric_rasters,
                                      c(arguments, list(n_trees = 15))))

  expect_equal(names(point$pd_stack), c("pd_time", "pdi_time", "pd_traits", "pdi_traits"))

  expect_equal(names(sampled$pd_stack),
               c("pd_time", "pdi_time",
                 "pd_traits_mean", "pd_traits_sd", "pd_traits_q2.5", "pd_traits_q97.5",
                 "pdi_traits_mean", "pdi_traits_sd", "pdi_traits_q2.5", "pdi_traits_q97.5"))

  # The time-based metrics do not involve the traits, so sampling must not move
  # them at all.
  expect_equal(raster::values(sampled$pd_stack[["pd_time"]]),
               raster::values(point$pd_stack[["pd_time"]]))
  expect_equal(raster::values(sampled$pd_stack[["pdi_time"]]),
               raster::values(point$pd_stack[["pdi_time"]]))

  expect_s3_class(sampled$trait_phylo, "multiPhylo")
  expect_length(sampled$trait_phylo, 15)
  expect_s3_class(point$trait_phylo, "phylo")

})

test_that("the summary layers are internally consistent and carry real spread", {

  skip_on_cran()

  fixture <- simulate_occurrences()

  template <- raster::raster(nrows = 10, ncols = 10)

  sampled <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = fixture$traits,
                                template.raster = template,
                                n_trees = 25))

  values <- raster::values(sampled$pd_stack)

  occupied <- !is.na(values[, "pd_traits_mean"])

  expect_gt(sum(occupied), 0)

  expect_true(all(values[occupied, "pd_traits_q2.5"] <= values[occupied, "pd_traits_mean"]))
  expect_true(all(values[occupied, "pd_traits_mean"] <= values[occupied, "pd_traits_q97.5"]))

  # If this were zero the uncertainty would not be reaching the metrics at all.
  expect_true(all(values[occupied, "pd_traits_sd"] > 0))

  # The same cells are NA in the summary layers as in the point estimate.
  point <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = fixture$traits,
                                template.raster = template))

  expect_equal(is.na(values[, "pd_traits_mean"]),
               is.na(raster::values(point$pd_stack[["pd_traits"]])))

})

test_that("probs controls the quantiles reported", {

  skip_on_cran()

  fixture <- simulate_occurrences()

  template <- raster::raster(nrows = 10, ncols = 10)

  sampled <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = fixture$traits,
                                template.raster = template,
                                n_trees = 15,
                                probs = c(0.25, 0.75)))

  expect_true(all(c("pd_traits_q25", "pd_traits_q75") %in% names(sampled$pd_stack)))

  values <- raster::values(sampled$pd_stack)
  occupied <- !is.na(values[, "pd_traits_mean"])

  expect_true(all(values[occupied, "pd_traits_q25"] <= values[occupied, "pd_traits_q75"]))

})

test_that("an already scaled phylogeny can be supplied instead of traits", {

  skip_on_cran()

  fixture <- simulate_occurrences()

  template <- raster::raster(nrows = 10, ncols = 10)

  internal <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = fixture$traits,
                                template.raster = template))

  supplied <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = NULL,
                                template.raster = template,
                                trait_phylo = internal$trait_phylo))

  expect_equal(raster::values(supplied$pd_stack[["pd_traits"]]),
               raster::values(internal$pd_stack[["pd_traits"]]))

  # No traits, so no rate-scaled phylogeny can be built.
  expect_null(supplied$rate_phylo)

  # A multiPhylo is summarised the same way n_trees would be.
  posterior <- suppressMessages(
    scale_branches_multidimensional_with_variation(fixture$tree, fixture$traits, n_trees = 12))

  from_sample <- suppressMessages(
    calculate_pd_metric_rasters(occurrences = fixture$occurrences,
                                phylogeny = fixture$tree,
                                traits = NULL,
                                template.raster = template,
                                trait_phylo = posterior))

  expect_true("pd_traits_mean" %in% names(from_sample$pd_stack))
  expect_length(from_sample$trait_phylo, 12)

})

test_that("conflicting or missing inputs are rejected", {

  fixture <- simulate_occurrences(n_tips = 15, n_cells = 15, n_records = 60)

  template <- raster::raster(nrows = 4, ncols = 4)

  expect_error(calculate_pd_metric_rasters(fixture$occurrences, fixture$tree,
                                           traits = NULL, template.raster = template),
               "Supply either traits")

  expect_error(calculate_pd_metric_rasters(fixture$occurrences, fixture$tree, fixture$traits,
                                           template, n_trees = 0),
               "n_trees must be a single positive number")

  expect_error(calculate_pd_metric_rasters(fixture$occurrences, fixture$tree, fixture$traits,
                                           template, probs = 0.5),
               "probs must be two probabilities")

  expect_error(calculate_pd_metric_rasters(fixture$occurrences, fixture$tree, fixture$traits,
                                           template, probs = c(-1, 2)),
               "probs must be two probabilities")

  # trait_phylo already fixes the trees, so n_trees would be ignored silently.
  scaled <- suppressMessages(
    scale_branches_multidimensional(fixture$tree, fixture$traits))

  expect_error(calculate_pd_metric_rasters(fixture$occurrences, fixture$tree, fixture$traits,
                                           template, trait_phylo = scaled, n_trees = 5),
               "n_trees applies only when")

  expect_error(calculate_pd_metric_rasters(fixture$occurrences, fixture$tree, fixture$traits,
                                           template, trait_phylo = "not a tree"),
               "must be a phylo object")

})
