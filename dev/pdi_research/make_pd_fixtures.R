# Freezes reference values from PhyloMeasures 2.1 so that the tests for the
# in-package PD and PDI code never need PhyloMeasures installed. Run from the
# package root, with PhyloMeasures available, only if the fixtures ever need
# regenerating:
#   remotes::install_url("https://cran.r-project.org/src/contrib/Archive/PhyloMeasures/PhyloMeasures_2.1.tar.gz")
#
# Cases 1 to 4 come straight from PhyloMeasures and can be regenerated at any
# time. The pipeline block is different: it is the output of the last version
# of calculate_pd_metric_rasters() that called PhyloMeasures (commit c5cf3df).
# Recomputing it with the current package would replace the independent
# reference with the implementation under test, so it is carried over from the
# existing fixture file, never recomputed here.

library(ape)

fixture_path <- "tests/testthat/fixtures/pd_reference.rds"

load_example <- function(name){

  env <- new.env()
  load(file.path("data", paste0(name, ".rda")), envir = env)
  get(name, envir = env)

}

example_tree <- load_example("example_tree")

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

# 5. end to end: the raster pipeline as computed through PhyloMeasures, taken
#    from the existing fixture and never recomputed (see the header)
if(!file.exists(fixture_path)){

  stop("The pipeline block can only come from the existing fixture file, which ",
       "is missing. It was produced by commit c5cf3df; check that out to rebuild it.")

}

pipeline <- readRDS(fixture_path)$pipeline

stopifnot(identical(names(pipeline), c("trait_tree", "trait_sample", "single", "several")))

saveRDS(list(ultrametric = ultrametric, ragged = ragged_case, tiny = tiny_case,
             ladder = ladder_case, pipeline = pipeline,
             source = "PhyloMeasures 2.1, null.model = 'uniform'"),
        fixture_path, version = 2)
