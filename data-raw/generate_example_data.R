# Generates the simulated example data shipped with the package.
# Run with source("data-raw/generate_example_data.R") from the package root.

library(ape)
library(phytools)

set.seed(20260820)

n_species <- 60

# A time-scaled phylogeny with binomial-style tip labels.
example_tree <- ape::rcoal(n = n_species)
example_tree$edge.length <- example_tree$edge.length / max(ape::node.depth.edgelength(example_tree))

genera <- c("Aphelia", "Brachyx", "Cnemodon", "Dolichor", "Eremita", "Fulgora")

example_tree$tip.label <- paste0(rep(genera, length.out = n_species), "_",
                                 sprintf("sp%02d", seq_len(n_species)))

# Three correlated traits evolving under Brownian motion, then put onto
# plausible measurement scales so that they differ in units, as real trait data
# do.  All are strictly positive.
correlations <- matrix(c(1.0, 0.6, 0.3,
                         0.6, 1.0, 0.4,
                         0.3, 0.4, 1.0), nrow = 3)

latent <- phytools::sim.corrs(tree = example_tree, vcv = correlations)

example_traits <- data.frame(
  species     = rownames(latent),
  body_mass_g = round(exp(3.2 + 0.8 * latent[, 1]), 2),
  wing_mm     = round(exp(4.1 + 0.3 * latent[, 2]), 2),
  bill_mm     = round(exp(2.3 + 0.25 * latent[, 3]), 2),
  row.names   = NULL,
  stringsAsFactors = FALSE
)

# Occurrences in tidy format over a 10 x 10 grid, with species occupying
# contiguous-ish sets of cells rather than being scattered uniformly.
n_cells <- 100

occurrence_rows <- do.call(rbind, lapply(seq_len(n_species), function(i) {

  centre <- sample(n_cells, 1)
  range_size <- sample(3:12, 1)

  cells <- unique(pmin(n_cells, pmax(1, centre + sample(-15:15, range_size, replace = TRUE))))

  data.frame(species = example_tree$tip.label[i],
             cell = cells,
             stringsAsFactors = FALSE)

}))

example_occurrences <- occurrence_rows[order(occurrence_rows$cell,
                                             occurrence_rows$species), ]
row.names(example_occurrences) <- NULL

stopifnot(setequal(example_traits$species, example_tree$tip.label),
          all(example_occurrences$species %in% example_tree$tip.label),
          !anyNA(example_traits))

usethis::use_data(example_tree, overwrite = TRUE)
usethis::use_data(example_traits, overwrite = TRUE)
usethis::use_data(example_occurrences, overwrite = TRUE)
