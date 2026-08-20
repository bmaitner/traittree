# scale_branches_multidimensional_with_variation() draws trees from the joint
# posterior of the ancestral states.  Helpers live in helper-branch_lengths.R
# and helper-reference-implementations.R.

test_that("observed tip values are carried through the sample unchanged", {

  sim <- simulate_traits(n_tips = 25, seed = 91)

  fit <- Rphylopars::phylopars(trait_data = sim$traits, tree = sim$tree, pheno_error = FALSE)

  observed <- as.matrix(sim$traits[match(sim$tree$tip.label, sim$traits$species), -1])

  for(i in 1:5){

    sampled <- sample_ancestral_states(tree = sim$tree,
                                       traits = sim$traits,
                                       posterior_mean = unname(fit$anc_recon),
                                       phylocov = fit$pars$phylocov,
                                       phenocov = fit$pars$phenocov,
                                       pheno_error = FALSE)

    # One row per node, tips first: uncertainty belongs to the unobserved nodes.
    expect_equal(nrow(sampled), 2 * length(sim$tree$tip.label) - 1)
    expect_equal(unname(sampled[seq_len(nrow(observed)), ]), unname(observed),
                 tolerance = 1e-8)

  }

})

test_that("draws recover the posterior mean and the per-node marginal variances", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 20, n_traits = 2, seed = 92)

  fit <- Rphylopars::phylopars(trait_data = sim$traits, tree = sim$tree, pheno_error = FALSE)

  n_draws <- 3000
  draws <- array(NA_real_, dim = c(n_draws, nrow(fit$anc_recon), ncol(fit$anc_recon)))

  set.seed(3)
  for(i in seq_len(n_draws)){

    draws[i, , ] <- sample_ancestral_states(tree = sim$tree,
                                            traits = sim$traits,
                                            posterior_mean = unname(fit$anc_recon),
                                            phylocov = fit$pars$phylocov,
                                            phenocov = fit$pars$phenocov,
                                            pheno_error = FALSE)

  }

  empirical_mean <- apply(draws, c(2, 3), mean)
  empirical_var <- apply(draws, c(2, 3), stats::var)

  # Monte Carlo error scales as sd/sqrt(n_draws); compare on that scale.
  standard_error <- sqrt(fit$anc_var/n_draws)

  expect_lt(max(abs(empirical_mean - fit$anc_recon)/pmax(standard_error, 1e-12)), 5)

  # Marginal variances are anc_var, which the fitted object reports directly.
  internal <- (length(sim$tree$tip.label) + 1):nrow(fit$anc_recon)
  expect_lt(max(abs(empirical_var[internal, ] - fit$anc_var[internal, ])/
                  fit$anc_var[internal, ]), 0.15)

})

test_that("draws reproduce the exact joint covariance between nodes", {

  skip_on_cran()

  # This is the property that marginal-only sampling gets wrong: branch lengths
  # depend on differences between adjacent nodes, so the covariance between
  # nodes has to be right, not just each node's own variance.

  sim <- simulate_traits(n_tips = 14, n_traits = 1, seed = 93)

  fit <- Rphylopars::phylopars(trait_data = sim$traits, tree = sim$tree, pheno_error = FALSE)

  tip_values <- stats::setNames(sim$traits$trait1, sim$traits$species)

  exact <- exact_joint_posterior(tree = sim$tree,
                                 x_tips = tip_values,
                                 R = as.numeric(fit$pars$phylocov))

  n_draws <- 4000
  internal <- (length(sim$tree$tip.label) + 1):nrow(fit$anc_recon)
  draws <- matrix(NA_real_, nrow = n_draws, ncol = length(internal))

  set.seed(4)
  for(i in seq_len(n_draws)){

    draws[i, ] <- sample_ancestral_states(tree = sim$tree,
                                          traits = sim$traits,
                                          posterior_mean = unname(fit$anc_recon),
                                          phylocov = fit$pars$phylocov,
                                          phenocov = fit$pars$phenocov,
                                          pheno_error = FALSE)[internal, 1]

  }

  expect_lt(max(abs(stats::cov(draws) - exact$cov))/max(abs(exact$cov)), 0.12)
  expect_lt(max(abs(colMeans(draws) - exact$mean))/max(sqrt(diag(exact$cov))), 0.15)

})

test_that("sampled branch lengths are not inflated the way marginal sampling would be", {

  skip_on_cran()

  # Drawing each node from its own marginal ignores the correlation between
  # neighbouring nodes and inflates branch lengths, worst on short branches.
  # The joint sampler should track the exact expectation instead.

  sim <- simulate_traits(n_tips = 30, n_traits = 1, seed = 94)

  fit <- Rphylopars::phylopars(trait_data = sim$traits, tree = sim$tree, pheno_error = FALSE)

  n_draws <- 400

  set.seed(5)
  joint <- suppressMessages(
    scale_branches_multidimensional_with_variation(tree = sim$tree,
                                                   traits = sim$traits,
                                                   n_trees = n_draws))

  joint_means <- colMeans(do.call(rbind, lapply(joint, function(x) x$edge.length)))

  # the same thing done from marginals only
  marginal <- matrix(NA_real_, nrow = n_draws, ncol = nrow(sim$tree$edge))
  n_tips <- length(sim$tree$tip.label)

  set.seed(6)
  for(i in seq_len(n_draws)){

    values <- stats::rnorm(nrow(fit$anc_recon), fit$anc_recon[, 1], sqrt(fit$anc_var[, 1]))
    values[seq_len(n_tips)] <- fit$anc_recon[seq_len(n_tips), 1]
    marginal[i, ] <- abs(values[sim$tree$edge[, 1]] - values[sim$tree$edge[, 2]])

  }

  marginal_means <- colMeans(marginal)

  shortest <- sim$tree$edge.length <= stats::quantile(sim$tree$edge.length, 0.25)

  # Marginal sampling should be visibly larger, especially on short branches.
  expect_gt(mean(marginal_means/joint_means), 1.1)
  expect_gt(mean((marginal_means/joint_means)[shortest]), 1.2)

})

test_that("n_trees controls the return type and the number of trees", {

  sim <- simulate_traits(n_tips = 20, seed = 95)

  one <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, sim$traits))

  expect_s3_class(one, "phylo")
  expect_length(one$edge.length, nrow(sim$tree$edge))
  expect_true(all(is.finite(one$edge.length)))
  expect_null(dim(one$edge.length))

  many <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, sim$traits, n_trees = 4))

  expect_s3_class(many, "multiPhylo")
  expect_length(many, 4)

  # Independent draws, so no two trees should coincide.
  lengths <- lapply(many, function(x) x$edge.length)
  expect_false(isTRUE(all.equal(lengths[[1]], lengths[[2]])))
  expect_false(isTRUE(all.equal(lengths[[3]], lengths[[4]])))

  expect_error(scale_branches_multidimensional_with_variation(sim$tree, sim$traits, n_trees = 0),
               "n_trees must be a single positive number")
  expect_error(scale_branches_multidimensional_with_variation(sim$tree, sim$traits, n_trees = "2"),
               "n_trees must be a single positive number")

})

test_that("rate = TRUE divides by the original branch lengths", {

  sim <- simulate_traits(n_tips = 20, seed = 96)

  set.seed(7)
  absolute <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, sim$traits, n_trees = 50))

  set.seed(7)
  rates <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, sim$traits, n_trees = 50, rate = TRUE))

  for(i in seq_along(absolute)){

    expect_equal(rates[[i]]$edge.length,
                 absolute[[i]]$edge.length/sim$tree$edge.length,
                 tolerance = 1e-10)

  }

})

test_that("missing trait values are sampled too, and observed ones are preserved", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, seed = 97)

  traits <- add_missing_values(sim$traits, proportion = 0.2)

  fit <- Rphylopars::phylopars(trait_data = traits, tree = sim$tree, pheno_error = FALSE)

  observed <- as.matrix(traits[match(sim$tree$tip.label, traits$species), -1])
  is_observed <- !is.na(observed)

  first <- sample_ancestral_states(tree = sim$tree, traits = traits,
                                   posterior_mean = unname(fit$anc_recon),
                                   phylocov = fit$pars$phylocov,
                                   phenocov = fit$pars$phenocov,
                                   pheno_error = FALSE)

  second <- sample_ancestral_states(tree = sim$tree, traits = traits,
                                    posterior_mean = unname(fit$anc_recon),
                                    phylocov = fit$pars$phylocov,
                                    phenocov = fit$pars$phenocov,
                                    pheno_error = FALSE)

  tip_rows <- seq_len(nrow(observed))

  # observed entries are fixed across draws...
  expect_equal(first[tip_rows, ][is_observed], observed[is_observed], tolerance = 1e-8)
  expect_equal(second[tip_rows, ][is_observed], observed[is_observed], tolerance = 1e-8)

  # ...and the imputed ones genuinely vary.
  expect_false(isTRUE(all.equal(first[tip_rows, ][!is_observed],
                                second[tip_rows, ][!is_observed])))

  scaled <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, traits, n_trees = 2))

  expect_length(scaled, 2)
  expect_true(all(is.finite(scaled[[1]]$edge.length)))

})

test_that("replicated observations run through the pheno_error path", {

  skip_on_cran()

  sim <- simulate_traits(n_tips = 30, seed = 98)

  set.seed(8)
  extra <- do.call(rbind, replicate(3, sim$traits[1:12, ], simplify = FALSE))
  extra[, -1] <- extra[, -1] +
    matrix(stats::rnorm(nrow(extra) * (ncol(extra) - 1), sd = 0.1), nrow = nrow(extra))

  replicated <- rbind(sim$traits, extra)

  expect_true(resolve_pheno_error(replicated))

  scaled <- suppressMessages(
    scale_branches_multidimensional_with_variation(sim$tree, replicated, n_trees = 3))

  expect_s3_class(scaled, "multiPhylo")
  expect_length(scaled, 3)

  for(tree_i in scaled){

    expect_length(tree_i$edge.length, nrow(sim$tree$edge))
    expect_true(all(is.finite(tree_i$edge.length)))

  }

})
