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
