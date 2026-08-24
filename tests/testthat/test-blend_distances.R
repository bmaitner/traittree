# blend_distances() and tune_blend_weight(): the alternative route to combining
# trait and phylogenetic information, for comparison against branch rescaling.

test_that("the endpoints reduce to the scaled inputs", {

  set.seed(1)
  tree <- ape::rtree(30)
  traits <- matrix(stats::rnorm(90), nrow = 30,
                   dimnames = list(tree$tip.label, c("t1", "t2", "t3")))

  phylo <- ape::cophenetic.phylo(tree)
  trait <- as.matrix(stats::dist(traits))
  species <- row.names(phylo)

  expect_equal(blend_distances(tree, traits, a = 1),
               (phylo / max(phylo))[species, species])

  expect_equal(blend_distances(tree, traits, a = 0),
               (trait / max(trait))[species, species])

})

test_that("the p-norm is applied as defined, and p = 2 differs from p = 1", {

  set.seed(2)
  tree <- ape::rtree(20)
  traits <- matrix(stats::rnorm(40), nrow = 20,
                   dimnames = list(tree$tip.label, c("t1", "t2")))

  phylo <- ape::cophenetic.phylo(tree); phylo <- phylo / max(phylo)
  trait <- as.matrix(stats::dist(traits)); trait <- trait / max(trait)
  species <- row.names(phylo)
  a <- 0.3

  expect_equal(blend_distances(tree, traits, a = a, p = 1),
               ((a * phylo) + ((1 - a) * trait))[species, species])

  expect_equal(blend_distances(tree, traits, a = a, p = 2),
               sqrt((a * phylo^2) + ((1 - a) * trait^2))[species, species])

  # the default is p = 2, and it is not the same thing as a linear mixture
  expect_false(isTRUE(all.equal(blend_distances(tree, traits, a = a, p = 1),
                                blend_distances(tree, traits, a = a, p = 2))))

})

test_that("the result is a valid distance matrix", {

  set.seed(3)
  tree <- ape::rtree(25)
  traits <- matrix(stats::rnorm(100), nrow = 25,
                   dimnames = list(tree$tip.label, paste0("t", 1:4)))

  blended <- blend_distances(tree, traits, a = 0.4)

  expect_equal(blended, t(blended))
  expect_equal(unname(diag(blended)), rep(0, nrow(blended)))
  expect_true(all(blended >= 0))

})

test_that("species are matched by name rather than by position", {

  set.seed(4)
  tree <- ape::rtree(15)
  traits <- matrix(stats::rnorm(30), nrow = 15,
                   dimnames = list(tree$tip.label, c("t1", "t2")))

  shuffled <- traits[sample(nrow(traits)), , drop = FALSE]

  expect_equal(blend_distances(tree, traits, a = 0.5),
               blend_distances(tree, shuffled, a = 0.5))

})

test_that("blend_distances rejects weights outside [0, 1] and too few shared species", {

  set.seed(5)
  tree <- ape::rtree(10)
  traits <- matrix(stats::rnorm(20), nrow = 10,
                   dimnames = list(tree$tip.label, c("t1", "t2")))

  expect_error(blend_distances(tree, traits, a = 1.5))
  expect_error(blend_distances(tree, traits, a = -0.1))

  unrelated <- matrix(stats::rnorm(20), nrow = 10,
                      dimnames = list(paste0("other", 1:10), c("t1", "t2")))
  expect_error(blend_distances(tree, unrelated, a = 0.5), "fewer than 3 species")

})

test_that("tune_blend_weight returns a weight on the grid and needs two traits", {

  set.seed(6)
  tree <- ape::rtree(40)
  traits <- matrix(stats::rnorm(120), nrow = 40,
                   dimnames = list(tree$tip.label, c("t1", "t2", "t3")))

  fit <- tune_blend_weight(tree, traits, a_grid = seq(0, 1, by = 0.25))

  expect_true(fit$a %in% fit$a_grid)
  expect_length(fit$mean_score, length(fit$a_grid))
  expect_equal(dim(fit$fold_scores), c(length(fit$a_grid), ncol(traits)))

  expect_error(tune_blend_weight(tree, traits[, 1, drop = FALSE]),
               "at least 2 traits")

})

test_that("tuning gives the phylogeny weight when it is informative", {

  # Traits simulated ON the tree are predictable from it, so the phylogeny should
  # earn weight and the chosen blend should beat traits alone.
  set.seed(7)
  tree <- ape::rtree(80)

  structured <- replicate(n = 4,
                          expr = ape::rTraitCont(phy = tree, model = "BM"))
  row.names(structured) <- tree$tip.label

  fit <- tune_blend_weight(tree, structured)

  expect_gt(fit$a, 0)
  expect_gt(max(fit$mean_score), fit$mean_score[fit$a_grid == 0])

})

test_that("the criterion is flat when nothing predicts the held-out trait", {

  # With traits that are independent of each other AND of the tree, no weight is
  # better than any other and the choice is arbitrary. Documented rather than
  # asserted away: a user whose score curve is flat should not read the returned
  # weight as meaningful.
  set.seed(8)
  tree <- ape::rtree(80)
  noise <- matrix(stats::rnorm(320), nrow = 80,
                  dimnames = list(tree$tip.label, paste0("t", 1:4)))

  fit <- tune_blend_weight(tree, noise)

  expect_lt(diff(range(fit$mean_score)), 0.15)

})

test_that("both fold criteria run and return comparable structures", {

  set.seed(8)
  tree <- ape::rtree(30)
  traits <- matrix(stats::rnorm(90), nrow = 30,
                   dimnames = list(tree$tip.label, c("t1", "t2", "t3")))

  mantel <- tune_blend_weight(tree, traits, a_grid = seq(0, 1, by = 0.25),
                              method = "mantel")
  unique <- tune_blend_weight(tree, traits, a_grid = seq(0, 1, by = 0.25),
                              method = "uniqueness")

  expect_equal(mantel$method, "mantel")
  expect_equal(unique$method, "uniqueness")
  expect_equal(dim(mantel$fold_scores), dim(unique$fold_scores))

})

test_that("p is validated: Inf is rejected and p < 1 warns", {

  set.seed(9)
  tree <- ape::rtree(15)
  traits <- matrix(stats::rnorm(30), nrow = 15,
                   dimnames = list(tree$tip.label, c("t1", "t2")))

  # the Chebyshev limit is a component-wise maximum in which `a` drops out, and
  # the p-norm expression evaluates the zero diagonal as 0^0, so it silently
  # returned a matrix of ones
  expect_error(blend_distances(tree, traits, a = 0.5, p = Inf), "finite")
  expect_error(blend_distances(tree, traits, a = 0.5, p = 0))
  expect_error(tune_blend_weight(tree, traits, p = Inf), "finite")

  # allowed, since it scored well in testing, but no longer a metric
  expect_warning(blend_distances(tree, traits, a = 0.5, p = 0.5),
                 "triangle inequality")

})

test_that("a_grid is validated, and a single candidate is allowed", {

  set.seed(10)
  tree <- ape::rtree(20)
  traits <- matrix(stats::rnorm(60), nrow = 20,
                   dimnames = list(tree$tip.label, c("t1", "t2", "t3")))

  expect_error(tune_blend_weight(tree, traits, a_grid = numeric(0)), "non-empty")
  expect_error(tune_blend_weight(tree, traits, a_grid = c(0.5, 1.5)), "between 0 and 1")
  expect_error(tune_blend_weight(tree, traits, a_grid = c(0.5, NA)), "between 0 and 1")

  # a one-point grid used to reach rowMeans() with a vector rather than a matrix
  fit <- tune_blend_weight(tree, traits, a_grid = 0.5)

  expect_equal(fit$a, 0.5)
  expect_equal(dim(fit$fold_scores), c(1L, ncol(traits)))

})

test_that("inputs with no variation stay usable instead of becoming NaN", {

  set.seed(11)
  tree <- ape::rtree(12)
  traits <- matrix(stats::rnorm(24), nrow = 12,
                   dimnames = list(tree$tip.label, c("t1", "t2")))

  constant <- matrix(1, nrow = 12, ncol = 2,
                     dimnames = list(tree$tip.label, c("t1", "t2")))
  flat_tree <- tree
  flat_tree$edge.length <- rep(0, nrow(tree$edge))

  # dividing by a zero maximum turned the whole blend into NaN, including at the
  # endpoints where the degenerate input carries no weight at all
  expect_true(all(is.finite(blend_distances(tree, constant, a = 1))))
  expect_true(all(is.finite(blend_distances(flat_tree, traits, a = 0))))

  # a constant held-out trait scores NA for its own fold, quietly, and the
  # remaining folds still decide the weight
  part_constant <- cbind(traits, t3 = 1)
  fit <- expect_silent(tune_blend_weight(tree, part_constant,
                                         a_grid = seq(0, 1, by = 0.25)))

  expect_true(all(is.na(fit$fold_scores[, "t3"])))
  expect_true(is.finite(fit$a))

})

test_that("a square trait table is not mistaken for a distance matrix", {

  set.seed(12)
  tree <- ape::rtree(4)

  # square and symmetric, but neither zero on the diagonal nor a distance matrix
  square_traits <- matrix(c(3, 1, 2, 4,
                            1, 5, 6, 2,
                            2, 6, 7, 1,
                            4, 2, 1, 8), nrow = 4, byrow = TRUE,
                          dimnames = list(tree$tip.label, paste0("t", 1:4)))

  expect_equal(blend_distances(tree, square_traits, a = 0),
               local({
                 trait <- as.matrix(stats::dist(square_traits))
                 (trait / max(trait))[tree$tip.label, tree$tip.label]
               }))

})

test_that("a distance matrix named on one margin only is accepted", {

  set.seed(13)
  tree <- ape::rtree(10)
  traits <- matrix(stats::rnorm(20), nrow = 10,
                   dimnames = list(tree$tip.label, c("t1", "t2")))

  phylo <- ape::cophenetic.phylo(tree)
  colnames(phylo) <- NULL

  expect_equal(blend_distances(phylo, traits, a = 0.5),
               blend_distances(tree, traits, a = 0.5))

  unnamed <- unname(ape::cophenetic.phylo(tree))
  expect_error(blend_distances(unnamed, traits, a = 0.5), "species names")

})

test_that("the returned plateau brackets the chosen weight", {

  set.seed(14)
  tree <- ape::rtree(60)
  structured <- replicate(n = 3,
                          expr = ape::rTraitCont(phy = tree, model = "BM"))
  row.names(structured) <- tree$tip.label

  fit <- tune_blend_weight(tree, structured, a_grid = seq(0, 1, by = 0.1))

  expect_length(fit$a_plateau, 2)
  expect_lte(fit$a_plateau[1], fit$a)
  expect_gte(fit$a_plateau[2], fit$a)
  expect_length(fit$standard_error, length(fit$a_grid))

})
