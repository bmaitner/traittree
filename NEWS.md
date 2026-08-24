# traittree 0.0.0.9000

## New features

* `blend_distances()` and `tune_blend_weight()` add the other route to combining
  trait and phylogenetic information: rather than rescaling branch lengths, they
  mix a phylogenetic distance matrix with a trait distance matrix, following the
  functional-phylogenetic distance of Cadotte, Albert and Walker (2013). Which
  route is preferable depends on the question. Blending returns a distance
  matrix, so phylogenetic diversity and evolutionary distinctiveness are not
  defined on its output, whereas `scale_branches_multidimensional()` returns a
  phylogeny; but blending can describe a species' distinctiveness better when
  trait coverage is partial. Which route wins turns out to depend on the
  objective rather than on the data: accumulating diversity over a set
  (complementarity, surrogacy, reserve selection) favoured the trait-scaled
  tree, while ranking individual species by distinctiveness favoured a blend.
  `?blend_distances` now says so under "Choosing between the two routes".

* `tune_blend_weight()` estimates the blending weight from the data. The weight
  has no default and no established estimator in the literature, where the
  convention is to compute the blend across the whole range and read the curve
  as a diagnostic. Each measured trait is held out in turn and the weight chosen
  to reproduce its distance structure from the remaining traits and the
  phylogeny, on the reasoning that held-out measured traits stand in for the
  trait axes that were never measured. Inspect `mean_score`: when nothing
  predicts the held-out trait the curve is flat and the returned weight is
  arbitrary rather than estimated.

* `tune_blend_weight()` also returns `a_plateau`, the range of weights whose
  mean score falls within one standard error of the best, and `standard_error`
  for every candidate. The score curve is usually flat and the returned weight
  is biased low, so the plateau is a fairer summary of what the folds can
  actually distinguish than the single best weight is.

* `p` is documented as fixed by design rather than left as an unexamined
  default. Scoring a grid over `a` and `p` together showed `p` to be weakly
  identified: performance moved by only 0.04 to 0.07 in correlation across `p`
  from 0.1 to 5, cross-validation recovered the best `p` in 9% to 55% of
  subsets, and tuning `p` alongside `a` beat a fixed `p = 2` in about a fifth
  of them. Values below 1 remain available but now warn, since they break the
  triangle inequality and so return a dissimilarity rather than a metric.

* Guidance on the trait-scaled tree as the blend's phylogenetic component is
  now scoped to where it held. It helped in simulation, but on real mammal data
  the dated tree was the equal or better component at every trait count tried:
  what a phylogenetic component contributes to a blend is its independence from
  the measured traits, which rescaling gives up.


## Bug fixes

* `blend_distances()` and `tune_blend_weight()` reject `p = Inf` instead of
  returning a matrix of ones. The p-norm expression evaluated the zero diagonal
  as `0^0`, and the limit is a component-wise maximum in which `a` has no
  effect, so it was never the calculation the user asked for.

* An input with no variation, such as a constant trait or a tree with zero
  branch lengths, no longer turns the whole blend into `NaN`. Scaling by the
  maximum is skipped when that maximum is zero, so the endpoint blends stay
  usable and a constant held-out trait simply scores `NA` for its own fold.

* `tune_blend_weight()` accepts a one-point `a_grid`, which previously errored
  in `rowMeans()`, and rejects an empty or out-of-range grid, which previously
  returned `numeric(0)` as the chosen weight.

* Species names on a distance matrix are checked before they are matched on.
  Both margins must name the same species in the same order: the same names in
  a different order used to pass, and matching by name then reordered the
  columns away from the rows they belong to, returning an asymmetric "distance"
  matrix with a non-zero diagonal and no complaint. A species named twice is
  refused as well, since subsetting by name silently kept the first row and
  dropped the other species from the blend.

* A square species-by-trait table is no longer mistaken for a distance matrix.
  Squareness and symmetry alone were enough to skip the conversion to distances;
  the diagonal must now be zero and the entries non-negative as well. A distance
  matrix carrying species names on one margin only is named from the other
  rather than failing when the species are matched.

* Documentation cross-references render as links. The package now sets
  `Roxygen: list(markdown = TRUE)`, so `[foo()]` references in the blending
  functions and the example datasets no longer appear as literal brackets.


## Breaking changes

* `scale_branches_multidimensional()` and
  `scale_branches_multidimensional_with_variation()` gain `scale_traits`, which
  **defaults to `TRUE`**. Each trait is now divided by its standard deviation
  across species before distances are taken, so branch lengths no longer depend
  on the units the traits happen to be measured in. Previously, re-expressing one
  trait in different units changed the results. Set `scale_traits = FALSE` for
  the old behaviour, which remains appropriate when the traits are already on a
  common scale and their relative spread is itself meaningful.

* `scale_branches_multidimensional_with_variation()` now samples from the
  **joint** posterior of the ancestral states. Previously each node was redrawn
  independently for every branch that touched it, so a returned tree was a
  collection of per-branch draws rather than a coherent realisation. Sampled
  branch lengths will differ from previous versions, and a given `set.seed()`
  will not reproduce older output.

* All functions that call `Rphylopars::phylopars()` now decide `pheno_error`
  from the data rather than inheriting it. `phylopars()` turns the within-species
  error term on whenever the argument is left unset, but with a single
  observation per species there is nothing to estimate it from. This does not
  change results for single-observation data, where `phylopars()` fits no such
  term either way.

* `richness_raster()`'s documented return value was describing the raster
  pipeline's output. It returns a single raster of species counts.

## New features

* `scale_branches_multidimensional_with_variation()` gains `n_trees`, returning a
  `multiPhylo` of that many draws (a `phylo` when `n_trees = 1`, as before), so a
  distribution of trees can be generated in one call.

* `calculate_pd_metric_rasters()` gains `n_trees`, `probs` and `trait_phylo`.
  With `n_trees > 1` the trait-based layers are replaced by their mean, standard
  deviation and two quantiles across the posterior sample, carrying
  reconstruction uncertainty through to the maps. `trait_phylo` accepts an
  already-scaled `phylo` or `multiPhylo` instead of scaling internally.

* `scale_branches_multidimensional()` and
  `scale_branches_multidimensional_with_variation()` gain `weights`, applied to
  the squared trait differences.

* `pheno_error` is exposed on all functions that reconstruct ancestral states,
  to override the choice made from the data.

* `calculate_pd_metric_rasters()` gains `verbose`, replacing an unconditional
  `print()` of the loop counter.

* New warnings for situations where a calculation's assumptions do not hold:
  traits whose spreads differ by more than tenfold when used in their original
  units; `percent = TRUE` with ancestral values near or below zero; and a
  within-species error term switched on by very sparse replication, which is
  more often a duplicated row than a real repeat measurement.

* Simulated example data: `example_tree`, `example_traits` and
  `example_occurrences`.

* Every exported function now has runnable examples, and there is a vignette,
  `vignette("traittree")`.

## Performance

Branch-length assignment was quadratic in tree size across the package, and is
now vectorised. Measured on simulated data:

* `scale_branches_multidimensional()`: the assignment step went from 11.9 s to
  under 0.01 s at 12,000 tips; results identical.
* `scale_branches_by_traits_rphylopars()`: 7.8x faster at 1,600 tips, and more
  with additional traits, since the loop was nested inside a per-trait loop.
* `scale_branches_multidimensional_with_variation()`: 228x faster at 1,600 tips.
* `richness_raster()`: 55x faster on a 300 by 300 grid (268 s to 4.8 s).
* `calculate_pd_metric_rasters()`: 34x faster with 3,600 cells, by grouping
  occurrences in one pass and querying `PhyloMeasures` in blocks rather than one
  cell at a time.
* `scale_branches_by_traits_fastAnc()` also had a quadratic lookup removed,
  though `phytools::fastAnc()` dominates its runtime, so end-to-end timings are
  unchanged.
* `tune_blend_weight()` is roughly 1.5x faster (400 tips, 5 traits, the default
  grid). It called `blend_distances()` at every point of the grid, so it
  re-coerced, re-matched and re-scaled both matrices thousands of times over;
  the combination step is now shared with `blend_distances()` and the
  preparation happens once. The chosen weight and the fold scores are
  unchanged.

## Other

* `geiger` removed from `Imports`; it was unused.
* Added a test suite (272 assertions).
* `vignette("traittree")` gains a section on blending, covering the weight, the
  plateau, why `p` is fixed, and which of the two routes suits which objective.
