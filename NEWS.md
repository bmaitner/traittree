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
  trait coverage is partial. A trait-scaled tree often makes a better
  phylogenetic component than a dated one, so the two functions compose.

* `tune_blend_weight()` estimates the blending weight from the data. The weight
  has no default and no established estimator in the literature, where the
  convention is to compute the blend across the whole range and read the curve
  as a diagnostic. Each measured trait is held out in turn and the weight chosen
  to reproduce its distance structure from the remaining traits and the
  phylogeny, on the reasoning that held-out measured traits stand in for the
  trait axes that were never measured. Inspect `mean_score`: when nothing
  predicts the held-out trait the curve is flat and the returned weight is
  arbitrary rather than estimated.


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

## Other

* `geiger` removed from `Imports`; it was unused.
* Added a test suite (272 assertions).
