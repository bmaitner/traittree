# Freezes reference values from PhyloMeasures 2.1 so that the tests for the
# in-package PD and PDI code never need PhyloMeasures installed. Run from the
# package root, with PhyloMeasures available, only if the fixtures ever need
# regenerating:
#   remotes::install_url("https://cran.r-project.org/src/contrib/Archive/PhyloMeasures/PhyloMeasures_2.1.tar.gz")
#
# The end-to-end block at the bottom was produced with the last version of
# calculate_pd_metric_rasters() that called PhyloMeasures (commit c5cf3df); it
# cannot be regenerated from later versions and is kept as-is.

suppressMessages(devtools::load_all("."))
library(ape)

set.seed(20261106)

community_matrix <- function(tree, sizes){
  n <- length(tree$tip.label)
  m <- t(vapply(sizes, function(k){ v <- numeric(n); v[sample.int(n, k)] <- 1; v }, numeric(n)))
  colnames(m) <- tree$tip.label
  m
}

reference_case <- function(tree, sizes){
  m <- community_matrix(tree, sizes)
  n <- length(tree$tip.label)
  list(tree = tree,
       communities = m,
       pd = PhyloMeasures::pd.query(tree, m, standardize = FALSE),
       pdi = PhyloMeasures::pd.query(tree, m, standardize = TRUE),
       moments = PhyloMeasures::pd.moments(tree, sample.sizes = 0:n))
}

# 1. the shipped example tree: ultrametric, binary
ultrametric <- reference_case(example_tree, c(0, 1, 2, 3, 5, 8, 13, 20, 33, 45, 59, 60, 2, 2, 7, 7))

# 2. a non-ultrametric tree with polytomies and a few zero-length branches
ragged <- ape::rtree(35)
ragged$edge.length[sample.int(nrow(ragged$edge), 6)] <- 0
ragged <- ape::di2multi(ragged, tol = 1e-12)
ragged$edge.length[sample.int(nrow(ragged$edge), 3)] <- 0
ragged_case <- reference_case(ragged, c(0, 1, 2, 3, 4, 6, 10, 17, 25, 34, 35, 5, 5))

# 3. a very small tree, where every sample size is an edge case
tiny <- ape::rtree(5)
tiny_case <- reference_case(tiny, c(0, 1, 2, 3, 4, 5, 2, 3))

# 4. a caterpillar, the worst case for the nested-pair enumeration
ladder <- ape::stree(20, type = "left")
ladder$edge.length <- stats::runif(nrow(ladder$edge), 0.1, 2)
ladder_case <- reference_case(ladder, c(2, 3, 5, 10, 19, 20, 4))

# 5. end to end: the raster pipeline as computed through PhyloMeasures
trait_tree <- scale_branches_multidimensional(example_tree, example_traits)
set.seed(99)
trait_sample <- scale_branches_multidimensional_with_variation(example_tree, example_traits, n_trees = 4)
template <- raster::raster(nrows = 10, ncols = 10)

single <- calculate_pd_metric_rasters(example_occurrences, example_tree, traits = NULL,
                                      template.raster = template, trait_phylo = trait_tree)
several <- calculate_pd_metric_rasters(example_occurrences, example_tree, traits = NULL,
                                       template.raster = template, trait_phylo = trait_sample)

pipeline <- list(trait_tree = trait_tree,
                 trait_sample = trait_sample,
                 single = raster::values(single$pd_stack),
                 several = raster::values(several$pd_stack))

saveRDS(list(ultrametric = ultrametric, ragged = ragged_case, tiny = tiny_case,
             ladder = ladder_case, pipeline = pipeline,
             source = "PhyloMeasures 2.1, null.model = 'uniform'"),
        "tests/testthat/fixtures/pd_reference.rds", version = 2)
