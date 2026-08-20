##########

#' Simulate Brownian motion at every node of a tree
#'
#' Internal helper.  Draws an unconditional realization of multivariate
#' Brownian motion over the whole tree, with the root fixed at zero.  Increments
#' along each edge are drawn from `MVN(0, edge_length * phylocov)`.
#'
#' @param tree A phylogeny with branch lengths.
#' @param phylocov Evolutionary rate matrix, traits by traits.
#' @return Matrix of node values, one row per node in `ape` node order.
#' @noRd
simulate_bm_nodes <- function(tree, phylocov){

  n_nodes <- max(tree$edge)
  n_edges <- nrow(tree$edge)

  #One MVN(0, phylocov) draw per edge, then scaled to that edge's length.
  standard_normals <- matrix(data = stats::rnorm(n_edges * ncol(phylocov)),
                             nrow = n_edges,
                             ncol = ncol(phylocov))

  increments <- (standard_normals %*% chol(phylocov)) * sqrt(tree$edge.length)

  values <- matrix(data = 0, nrow = n_nodes, ncol = ncol(phylocov))

  #Cladewise order visits parents before children, so each node's value is
  #available by the time its descendants are reached.
  for(edge in ape::reorder.phylo(tree, order = "cladewise", index.only = TRUE)){

    values[tree$edge[edge, 2], ] <- values[tree$edge[edge, 1], ] + increments[edge, ]

  }

  values

}

#' Posterior mean trait values at every node
#'
#' Internal helper.  Returns the conditional expectation of every node's traits
#' given the supplied observations, holding the evolutionary parameters fixed.
#'
#' Complete data with one observation per species take a fast path through
#' `Rphylopars::anc.recon()`.  Missing values and replicate observations are
#' handled by `Rphylopars::phylopars()` with the covariance parameters fixed,
#' which skips re-estimation and only performs the reconstruction.
#'
#' @param tree A phylogeny with branch lengths.
#' @param traits A data frame whose first column is `species`.
#' @param phylocov Evolutionary rate matrix.
#' @param phenocov Within-species covariance matrix, or NULL.
#' @param pheno_error Whether a within-species error term is in the model.
#' @return Matrix of posterior means, one row per node in `ape` node order.
#' @noRd
reconstruct_all_nodes <- function(tree, traits, phylocov, phenocov, pheno_error){

  trait_columns <- setdiff(colnames(traits), "species")

  simple <- !pheno_error &&
    !anyNA(traits[, trait_columns, drop = FALSE]) &&
    !anyDuplicated(as.character(traits$species))

  if(simple){

    tip_values <- as.matrix(traits[match(tree$tip.label, traits$species),
                                   trait_columns, drop = FALSE])

    #anc.recon() looks its tips up by name.
    row.names(tip_values) <- tree$tip.label

    internal_values <- as.matrix(Rphylopars::anc.recon(trait_data = tip_values,
                                                       tree = tree))

    output <- rbind(unname(tip_values), unname(internal_values))

    colnames(output) <- trait_columns

    return(output)

  }

  arguments <- list(trait_data = traits,
                    tree = tree,
                    pheno_error = pheno_error,
                    phylocov_fixed = phylocov,
                    skip_optim = TRUE)

  if(pheno_error && !is.null(phenocov)){

    arguments$phenocov_fixed <- phenocov

  }

  unname(do.call(Rphylopars::phylopars, arguments)$anc_recon)

}

#' Draw one joint sample of ancestral (and missing tip) trait values
#'
#' Internal helper.  Returns a single draw from the joint posterior of every
#' node's traits given the observed data.
#'
#' Sampling is done by conditioning by kriging: an unconditional Brownian motion
#' realization is generated over the tree, reconstructed from its own simulated
#' observations, and the resulting residual is added to the posterior mean of the
#' real data. The result is an exact draw from the joint posterior, and costs one
#' reconstruction per sample rather than forming the covariance matrix over all
#' nodes, which would be quadratic in tree size.
#'
#' Sampling each node independently from its marginal variance would be wrong
#' here. Neighbouring nodes are strongly correlated, and branch lengths depend on
#' the difference between adjacent nodes, so ignoring that correlation inflates
#' branch lengths - severely for short branches.
#'
#' @param tree A phylogeny with branch lengths.
#' @param traits A data frame whose first column is `species`.
#' @param posterior_mean Posterior mean at every node, from the fitted model.
#' @param phylocov Evolutionary rate matrix.
#' @param phenocov Within-species covariance matrix, or NULL.
#' @param pheno_error Whether a within-species error term is in the model.
#' @return Matrix of sampled values, one row per node in `ape` node order.
#' @noRd
sample_ancestral_states <- function(tree,
                                    traits,
                                    posterior_mean,
                                    phylocov,
                                    phenocov,
                                    pheno_error){

  trait_columns <- setdiff(colnames(traits), "species")

  simulated_nodes <- simulate_bm_nodes(tree = tree, phylocov = phylocov)

  #Give the simulated data the same shape as the real data: the same species in
  #the same rows, the same pattern of missing values, and, when a within-species
  #error term is in the model, the same observation-level noise.
  simulated_traits <- traits

  tip_rows <- match(as.character(traits$species), tree$tip.label)

  simulated_values <- simulated_nodes[tip_rows, , drop = FALSE]

  if(pheno_error && !is.null(phenocov)){

    observation_noise <- matrix(data = stats::rnorm(nrow(simulated_values) * ncol(phenocov)),
                                nrow = nrow(simulated_values),
                                ncol = ncol(phenocov)) %*% chol(phenocov)

    simulated_values <- simulated_values + observation_noise

  }

  simulated_values[is.na(as.matrix(traits[, trait_columns, drop = FALSE]))] <- NA

  simulated_traits[, trait_columns] <- simulated_values

  simulated_mean <- reconstruct_all_nodes(tree = tree,
                                          traits = simulated_traits,
                                          phylocov = phylocov,
                                          phenocov = phenocov,
                                          pheno_error = pheno_error)

  posterior_mean + (simulated_nodes - simulated_mean)

}
