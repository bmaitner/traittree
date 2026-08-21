# traittree

Re-scale phylogenies from units of time into units of phenotypic change.

A time-scaled phylogeny measures branches in millions of years. For many
questions the more relevant currency is how much the organisms actually changed:
two lineages separated by the same span of time may have diverged a great deal
or hardly at all. `traittree` reconstructs ancestral trait values and rewrites
each branch length as the amount of trait change along it. Additional functions
propagate uncertainty into downstream products (e.g., phylogenetic diversity 
rasters).

## Installation

```r
# install.packages("remotes")
remotes::install_github("bmaitner/traittree")
```

## What it does

```r
library(traittree)

# Branch lengths become multivariate trait change rather than elapsed time
trait_tree <- scale_branches_multidimensional(tree = example_tree,
                                              traits = example_traits)

plot(example_tree$edge.length, trait_tree$edge.length,
     xlab = "time", ylab = "trait change")
```

Ancestral states are estimated, not observed, and that uncertainty can be
carried through rather than ignored. Each tree below is a draw from the joint
posterior of the ancestral states:

```r
trees <- scale_branches_multidimensional_with_variation(tree = example_tree,
                                                        traits = example_traits,
                                                        n_trees = 100)
```

Those trees can be pushed all the way through to the spatial metrics, so that
phylogenetic diversity comes with an interval rather than a single number:

```r
template <- raster::raster(nrows = 10, ncols = 10)

result <- calculate_pd_metric_rasters(occurrences = example_occurrences,
                                      phylogeny = example_tree,
                                      traits = example_traits,
                                      template.raster = template,
                                      n_trees = 100)

names(result$pd_stack)
#> "pd_time" "pdi_time" "pd_traits_mean" "pd_traits_sd" "pd_traits_q2.5" ...
```

## Functions

| Function | Purpose |
|---|---|
| `scale_branches_multidimensional()` | Branch lengths as multivariate trait distance |
| `scale_branches_multidimensional_with_variation()` | The same, drawn from the posterior of the ancestral states |
| `scale_branches_by_traits_rphylopars()` | One tree per trait |
| `scale_branches_by_traits_fastAnc()` | A single trait, via `phytools::fastAnc()` |
| `data_matching()` | Prune a tree and an occurrence table to their shared species |
| `richness_raster()` | Species richness per raster cell |
| `calculate_pd_metric_rasters()` | PD and PDI rasters from time- and trait-scaled trees |

## Two things worth knowing

**Traits are scaled by default.** Branch lengths are Euclidean distances across
traits, so a trait measured in grams would otherwise swamp one measured in
metres, and switching units would change the answer. Each trait is divided by
its standard deviation across species first. Set `scale_traits = FALSE` if your
traits are already on a common scale and their relative spread is itself
meaningful, as with log-transformed measurements.

**Within-species error is decided from the data.** `phylopars()` estimates a
within-species error term whenever it is left to its own devices, but with one
observation per species there is nothing to estimate it from. The package turns
it on only when at least one species has repeat observations of a trait, and
`pheno_error` overrides that either way.

## Vignette

```r
vignette("traittree", package = "traittree")
```

## License

MIT
