# Shared fixtures for the scale_branches_multidimensional() tests.

# The original edge-at-a-time loop, kept as a reference implementation so the
# vectorized replacement can be checked against it for both correctness and
# speed.  See issue #1.

reference_branch_lengths <- function(tree, anc_recon, rate = FALSE){

  row.names(anc_recon)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)

  output_branches <- numeric(length(tree$edge.length))

  for(i in 1:length(tree$edge.length)){

    node_1 <- tree$edge[i,][1]
    node_2 <- tree$edge[i,][2]

    value_1 <- anc_recon[which(row.names(anc_recon) == node_1),]
    value_2 <- anc_recon[which(row.names(anc_recon) == node_2),]

    bl <- stats::dist(rbind(value_1, value_2))[1]

    if(rate){

      bl <- bl/tree$edge.length[i]

    }

    output_branches[i] <- bl

  }

  output_branches

}

# Simulate a small tree with several correlated traits.

simulate_traits <- function(n_tips = 30, n_traits = 4, seed = 1){

  set.seed(seed)

  tree <- ape::rcoal(n = n_tips)
  tree$tip.label <- paste0("sp", 1:n_tips)

  # Correlated Brownian motion: a random covariance matrix shared across traits.
  loadings <- matrix(stats::rnorm(n_traits * n_traits), nrow = n_traits)
  sigma <- crossprod(loadings) + diag(n_traits)

  trait_values <- phytools::sim.corrs(tree = tree, vcv = sigma)

  # sim.corrs returns a bare named vector for a single trait.
  if(is.null(dim(trait_values))){

    trait_values <- matrix(trait_values,
                           ncol = 1,
                           dimnames = list(names(trait_values), NULL))

  }

  traits <- data.frame(species = row.names(trait_values),
                       trait_values,
                       row.names = NULL,
                       stringsAsFactors = FALSE)

  names(traits)[-1] <- paste0("trait", 1:n_traits)

  list(tree = tree, traits = traits)

}

# Punch NAs into every trait column, leaving each species some observed traits.

add_missing_values <- function(traits, proportion = 0.3, seed = 3){

  set.seed(seed)

  for(j in 2:ncol(traits)){

    traits[sample(x = nrow(traits), size = floor(proportion * nrow(traits))), j] <- NA

  }

  traits

}

# Median elapsed seconds for one evaluation of expr.  Repeats until at least
# `min_time` of work has accumulated, so sub-millisecond operations can be timed
# without the clock resolution dominating.

time_per_call <- function(expr, min_time = 0.2){

  expression <- substitute(expr)
  env <- parent.frame()

  repetitions <- 1

  repeat {

    elapsed <- system.time(for(i in seq_len(repetitions)) eval(expression, env))[["elapsed"]]

    if(elapsed >= min_time || repetitions > 1e5){

      return(elapsed/repetitions)

    }

    repetitions <- repetitions * 4

  }

}
