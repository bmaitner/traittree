# traittree

<!-- badges: start -->
[![R-CMD-check](https://github.com/bmaitner/traittree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bmaitner/traittree/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

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
# install.packages("pak")
pak::pak("bmaitner/traittree")
```

`remotes::install_github("bmaitner/traittree")` works equally well.

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

## Phylogenetic diversity, standardized exactly

The diversity calculations behind those maps are available on their own, for
any set of communities and any phylogeny:

```r
diversity <- phylo_diversity(tree = example_tree,
                             communities = example_occurrences,
                             rarefy_to = 4)

null <- phylo_diversity_moments(tree = example_tree, richness = 1:60)
```

`phylo_diversity()` returns Faith's PD and the phylogenetic diversity index
(PDI), which standardizes PD against the mean and standard deviation of every
set of the same richness drawn from the tree. Those moments are computed
exactly, from the formulas of Nipperess and Matsen (2013), rather than
estimated by shuffling, so the index carries no sampling error.

These replace `pd.query()` and `pd.moments()` from PhyloMeasures (Tsirogiannis
and Sandel 2016), which was archived from CRAN in 2022. Under its uniform null
model they return the same values, checked against values frozen from
PhyloMeasures to better than 1e-9. Beyond parity:

- **A rarefied index.** The standardized index still depends on richness:
  under the same degree of clustering, richer communities score more extremely
  (Sandel 2018). `rarefy_to` rarefies every community to a common richness. The
  average over rarefactions that Sandel estimates by replication is computed
  here exactly, in closed form.
- **Many trees at once.** Nothing in the null model depends on branch lengths,
  so trees that share a topology, such as a posterior sample from
  `scale_branches_multidimensional_with_variation()`, share all of the work.
  Across a sample of 25 to 50 trees this is roughly ten times faster than
  querying each tree in turn; for a single tree the speed is about the same.
- **Numerics.** Against a pairwise-covariance reference on a 1,500-tip tree, the
  null standard deviation is accurate to about one part in 1e11, one to two
  orders of magnitude closer than PhyloMeasures, and it is exactly zero where
  the null has no spread.
- **Plain inputs.** Communities as a presence/absence matrix or as a table of
  occurrences; no C++ to compile.

Two things PhyloMeasures offered are not here: its abundance-weighted null
models (`frequency.by.richness`, `sequential`), and Monte Carlo p-values. The
index itself assumes the null is roughly normal, which it is not at very low
richness; treat values near the conventional cutoffs with that in mind.

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
| `phylo_diversity()` | PD and its standardized index (PDI) for any set of communities, with exact null moments and optional rarefaction |
| `phylo_diversity_moments()` | The null mean and standard deviation of PD at a given richness |
| `blend_distances()` | Mix phylogenetic and trait distances into one distance matrix |
| `tune_blend_weight()` | Estimate the blending weight from the data |

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
