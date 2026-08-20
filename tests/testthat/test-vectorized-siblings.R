# The remaining loop-based branch assignments, vectorized.  Each is checked
# against the original implementation in helper-reference-implementations.R.

test_that("fastAnc branch lengths match the reference loop", {

  sim <- simulate_traits(n_tips = 60, n_traits = 1, seed = 31)

  # This function takes a named vector rather than a data frame.
  trait_vector <- sim$traits$trait1
  names(trait_vector) <- sim$traits$species

  for(percent in c(FALSE, TRUE)){

    # simulated traits are centred on zero, so percent = TRUE warns about the
    # denominator; that warning has its own tests.
    scaled <- suppressWarnings(
      scale_branches_by_traits_fastAnc(tree = sim$tree,
                                       traits = trait_vector,
                                       percent = percent))

    expected <- reference_scale_branches_by_traits_fastAnc(tree = sim$tree,
                                                           traits = trait_vector,
                                                           percent = percent)

    expect_equal(as.numeric(scaled$edge.length),
                 as.numeric(expected$edge.length),
                 tolerance = 1e-10)

    expect_s3_class(scaled, "phylo")
    expect_null(dim(scaled$edge.length))

  }

})

test_that("fastAnc still drops species missing from the trait data", {

  sim <- simulate_traits(n_tips = 60, n_traits = 1, seed = 31)

  trait_vector <- sim$traits$trait1
  names(trait_vector) <- sim$traits$species

  dropped <- trait_vector[1:50]

  scaled <- scale_branches_by_traits_fastAnc(tree = sim$tree, traits = dropped)

  expect_length(scaled$tip.label, 50)
  expect_setequal(scaled$tip.label, names(dropped))
  expect_length(scaled$edge.length, nrow(scaled$edge))

})

test_that("rphylopars branch lengths match the reference loop, trait by trait", {

  sim <- simulate_traits(n_tips = 40, seed = 32)

  for(percent in c(FALSE, TRUE)){

    scaled <- suppressWarnings(suppressMessages(
      scale_branches_by_traits_rphylopars(tree = sim$tree,
                                          traits = sim$traits,
                                          percent = percent)))

    expected <- reference_scale_branches_by_traits_rphylopars(tree = sim$tree,
                                                              traits = sim$traits,
                                                              percent = percent)

    expect_s3_class(scaled, "multiPhylo")
    expect_length(scaled, ncol(sim$traits) - 1)
    expect_equal(names(scaled), names(expected))

    for(x in seq_along(scaled)){

      expect_equal(as.numeric(scaled[[x]]$edge.length),
                   as.numeric(expected[[x]]$edge.length),
                   tolerance = 1e-10)

    }

    expect_null(dim(scaled[[1]]$edge.length))

  }

})

test_that("rphylopars variant resolves and forwards pheno_error", {

  sim <- simulate_traits(n_tips = 20, seed = 33)

  replicated <- rbind(sim$traits, sim$traits[1:3, ])

  captured <- NULL

  fake_phylopars <- function(trait_data, tree, pheno_error, ...){

    captured <<- pheno_error

    n_tips <- length(tree$tip.label)
    anc_recon <- matrix(seq_len((2 * n_tips - 1) * (ncol(trait_data) - 1)),
                        nrow = 2 * n_tips - 1)
    colnames(anc_recon) <- setdiff(colnames(trait_data), "species")
    row.names(anc_recon) <- c(tree$tip.label, as.character((n_tips + 1):(2 * n_tips - 1)))

    list(anc_recon = anc_recon)

  }

  local_mocked_bindings(phylopars = fake_phylopars, .package = "traittree")

  suppressMessages(scale_branches_by_traits_rphylopars(sim$tree, sim$traits))
  expect_false(captured)

  suppressMessages(scale_branches_by_traits_rphylopars(sim$tree, replicated))
  expect_true(captured)

  suppressMessages(scale_branches_by_traits_rphylopars(sim$tree, sim$traits, pheno_error = TRUE))
  expect_true(captured)

})

test_that("richness_raster matches the reference implementation", {


  template <- raster::raster(nrows = 20, ncols = 20)

  set.seed(35)
  occurrences <- data.frame(species = paste0("sp", sample(50, 400, replace = TRUE)),
                            cell = sample(400, 400, replace = TRUE),
                            stringsAsFactors = FALSE)

  # identical, not merely equal: counts should stay integers, as before.
  expect_identical(raster::values(richness_raster(template.raster = template,
                                                  occurrences = occurrences)),
                   raster::values(reference_richness_raster(template.raster = template,
                                                            occurrences = occurrences)))

})

test_that("richness_raster counts distinct species and leaves empty cells NA", {


  template <- raster::raster(nrows = 5, ncols = 5)

  # Cell 1: two distinct species, one of them recorded twice.
  # Cell 3: a single species.
  occurrences <- data.frame(species = c("a", "b", "a", "c"),
                            cell = c(1, 1, 1, 3),
                            stringsAsFactors = FALSE)

  values <- raster::values(richness_raster(template.raster = template,
                                           occurrences = occurrences))

  expect_equal(values[1], 2)
  expect_equal(values[3], 1)
  expect_true(all(is.na(values[c(2, 4:25)])))

})
