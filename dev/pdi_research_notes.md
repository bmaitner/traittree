# Notes: computing PD and PDI without PhyloMeasures

Research notes from November 2026, ahead of replacing the two
`PhyloMeasures::pd.query()` calls in `calculate_pd_metric_rasters()` so that the
package can go to CRAN. This folder is build-ignored.

## What we need to reproduce

`pd.query(tree, matrix, standardize = FALSE)` and `standardize = TRUE`, under
the default `null.model = "uniform"`. Read from the PhyloMeasures 2.1 source and
confirmed numerically:

* **PD** is the total branch length of the minimal subtree connecting the
  sample (unrooted PD, the Steiner tree). The path from that subtree up to the
  root is *not* included. Checked on a deep clade: 0.181 without the root path,
  1.143 with it; `pd.query` returns 0.181.
* **PDI** is `(PD - mean) / sd`, with the mean and sd taken over all samples of
  exactly the same richness drawn uniformly without replacement from the tips of
  the tree. This is eq. 1 of Sandel (2018). Negative means clustered.
* Edge cases in the C++: an empty sample returns 0; richness 1 gives PD 0, mean
  0, sd 0; when the sd is zero (richness 1, or richness equal to the pool) the
  standardized value returned is `PD - mean`, undivided, rather than NaN.
* The pool is every tip of the tree passed in, so pruning the tree changes the
  null. `calculate_pd_metric_rasters()` prunes to species with occurrences
  before querying; that behaviour has to be kept for results to match.

## The mathematics

Let the tree have `s` tips, the sample have `r`, edge `e` have length `w_e` and
`s_e` tips below it, and `h(x) = choose(x, r) / choose(s, r)` be the probability
that all `r` draws fall inside a given set of `x` tips (zero for `x < r`).

An edge is in the sampled subtree unless every draw is below it or every draw is
elsewhere, so

    P(e in) = 1 - h(s_e) - h(s - s_e)
    E[PD]   = sum_e w_e P(e in)

For the variance, `E[PD^2] = sum_e sum_f w_e w_f P(e and f both in)`, with

    clades disjoint:        1 - h(s - s_e) - h(s - s_f) + h(s - s_e - s_f)
    e inside f (s_e <= s_f): 1 - h(s_f)     - h(s - s_e) + h(s_f - s_e)
    e = f:                   P(e in)

and `Var[PD] = E[PD^2] - E[PD]^2`.

These are exactly Nipperess and Matsen (2013), eq. 5 for the mean and the
Appendix for the unrooted variance (their `q_k(R)`, the probability nothing in
`R` is drawn, is `h(s - |R|)`). The same structure is visible in PhyloMeasures'
`compute_variance_uniform_distribution_slow()`.

**Verified.** A brute-force evaluation of the formulas above over all edge pairs
reproduces `PhyloMeasures::pd.moments()` on `example_tree` for r = 2, 3, 8, 20,
45, 59, 60: means agree to 1e-14 and sds to 1e-7. The whole of the sd gap is at
r = s, where the true variance is zero and the subtraction cancels, so clamp the
variance at zero. Monte Carlo (20,000 draws, r = 8) agrees as well. Script:
`dev/pdi_research/verify_moments.R`; worth turning into a test.

## Cost, and what the fast algorithm does

Summing over all edge pairs is quadratic in tree size. The route PhyloMeasures
takes (Tsirogiannis, Sandel and Kalvisa 2014, as described in Tsirogiannis and
Sandel 2016, PLoS ONE) rests on two observations:

1. The disjoint-clade probability depends on the two edges only through their
   clade sizes, so edges can be pooled by clade size and the sum taken over
   pairs of size classes.
2. Treating every pair as disjoint is wrong only for ancestor-descendant pairs,
   which can be enumerated and corrected. There are as many of those as the sum
   of node depths.

The stated worst case is O(Sackin index), the sum of tip depths: O(n log n) for
balanced trees, O(n^2) for a caterpillar.

Something PhyloMeasures cannot exploit and we can: **the probabilities depend
only on topology and richness, never on branch lengths.** Every tree drawn by
`scale_branches_multidimensional_with_variation()` has the same topology. So for
a posterior sample the pair structure and probabilities are computed once per
distinct richness, and each tree costs one weighted sum. With `W` the edges by
trees matrix of branch lengths, the mean is `crossprod(p, W)` and the variance is
a quadratic form in each column.

## Why not picante

`picante::pd(include.root = FALSE)` matches `pd.query` exactly but is a pure R
loop: 69 times slower at 200 species by 1,600 cells.

`picante::expected.pd()` and `variance.pd()` are not the same quantity. Their
documentation states binomial sampling, each tip included independently with a
fixed probability, so richness is random and `n` is only its expectation; and
they are for rooted PD. On `example_tree` the means differ from PhyloMeasures by
+0.112, +0.007, -0.022 and -0.006 at r = 3, 8, 20, 45, and the sds are 25 to 50%
larger. For standardizing a community whose richness is known, fixed-size
sampling is the appropriate null, so these are the wrong tool rather than an
alternative convention.

## Sandel (2018): the index is still richness-dependent

Read in full. The standardization guarantees one thing only: communities
assembled at random score zero at any richness. It does not make non-zero values
comparable across richness. Holding the strength of an environmental filter
fixed, PDI became 3 to 4 times more negative as richness rose from 10 to 320
(rNRI about 8 times; rNTI 3 to 4 times), on 1,000-tip Yule trees.

His remedy is to rarefy every community to a common richness, the smallest in
the set, and compute the index on the rarefied communities. That removed the
dependence, slightly over-correcting. The cost is power: the rank correlation
between filter strength and PDI fell from 0.76 to 0.27 when rarefying 320
species to 10, and to 0.59 when rarefying to 80. Averaging over 1,000 replicate
rarefactions did not move the median but reduced the spread. His advice is to
prefer random error to a known directional bias, and to avoid extreme
rarefaction.

This matters here because the rasters compare cells of very different richness,
which is the situation the paper warns about.

**The replicated rarefaction has an exact closed form.** For a fixed target
richness `k` the index is an affine function of PD, so its average over
rarefactions is

    ( E_k[PD of the community's own subtree] - mu_k(T) ) / sigma_k(T)

and the expectation in the numerator is the same mean formula applied to the
subtree spanned by the community, with `s` replaced by the community's richness.
Checked against 4,000 replicate rarefactions to k = 10 at richness 15, 30 and
45: differences of -0.005, -0.013 and -0.008, each inside one Monte Carlo
standard error (`dev/pdi_research/verify_rarefied.R`). So Sandel's "RR" estimator can be had exactly, with no
replication, from machinery we need anyway. The per-edge counts it requires are
the ones the PD calculation produces.

Implemented as `rarefy_to` in `calculate_pd_metric_rasters()`; the tests check
it against the average over every possible rarefaction, which it matches to
rounding error.

## Other things learned

* Tsirogiannis and Sandel (2014, Algorithms for Molecular Biology) make the
  case that the null distribution of MPD is skewed, so reading a z-score against
  +/-1.96 gives biased p-values. The same caution applies to PD at low richness.
  PhyloMeasures gets PD p-values by Monte Carlo, not analytically.
* PhyloMeasures offers two further null models, `frequency.by.richness` and
  `sequential`, which weight species by abundance. We use neither.
* `R CMD check --as-cran` on the current package reports only NOTEs: the
  development version number, `PhyloMeasures` not in a mainstream repository,
  the Title not in title case, and the Description opening with "This package".

## Where it ended up

The internals are in `R/pd_internals.R`; `phylo_diversity()` and
`phylo_diversity_moments()` in `R/phylo_diversity.R` are the exported interface;
`calculate_pd_metric_rasters()` calls the internals directly.

## Deferred: abundance-weighted null models

PhyloMeasures offered two abundance-weighted nulls that the package does not
have. Deferred in November 2026 as not needed for the paper. What was learned
from the PhyloMeasures source, so that this can be picked up without re-reading
it:

**`sequential`.** Species drawn one at a time without replacement, each with
probability proportional to its abundance: what `sample(x, r, prob = w)` does.
PhyloMeasures has no formula for it and estimates the moments by Monte Carlo
(`reps` draws, default 1000). The same approach here would be cheap, since
`pd_of_communities()` handles thousands of replicate communities at once.

**`frequency.by.richness`.** Each species included independently with
probability `p_i = a_i / sum(a)`, the draw kept only if it has exactly `r`
species. Equivalent to `P(S) proportional to prod_{i in S} p_i / (1 - p_i)`
over sets of size `r`. PhyloMeasures computes the moments exactly with
generating functions (Tsirogiannis and Sandel 2016, RECOMB, "Fast phylogenetic
biodiversity computations under a non-uniform random distribution"; paywalled,
not read). From the source: each clade carries the polynomial
`prod_{tips i} ((1 - p_i) + p_i x)`, whose coefficient `k` is the probability
that `k` of its tips are selected; polynomials multiply up the tree, truncated
at degree `r`. For the mean, each edge contributes its length times the
coefficients `1..r-1` of its clade polynomial, and the cases where all `r`
selected tips lie below one edge (edge not on the unrooted subtree) are
subtracted using the probability that nothing is selected elsewhere. For the
variance, pairwise products of child polynomials across children, cross terms
of an edge with its descendants, and squared self terms. Everything is finally
divided by the coefficient `r` of the root polynomial. Two complications: the
coefficients underflow doubles on large trees (PhyloMeasures carries scaled
numbers; work in log scale here), and the pairwise variance terms are fiddly.
A day or two of careful work.

**Before starting either:** freeze reference values from PhyloMeasures with
`pd.query(..., null.model = ..., abundance.weights = ...)` and
`pd.moments(...)`, as `make_pd_fixtures.R` does for the uniform null. That
needs PhyloMeasures installed, which it is on this machine; it is the thing
that makes the implementation testable.

**Design decisions to make:** what the raster pipeline would use as
abundances (each species' number of occupied cells is the natural choice);
and that `rarefy_to` should stay restricted to the uniform null, since
rarefying a community is uniform subsampling and no reference defines its
combination with a weighted null.

## Licence

PhyloMeasures is GPL-3 and traittree is MIT, so its code cannot be translated
into this package. The formulas are published mathematics (Nipperess and Matsen
2013; Sandel 2018, eq. 1) and the implementation should be written from those,
in R, with its own structure. The C++ was read to settle conventions, not as a
template. If matching its structure ever seems necessary, ask the authors.

## Not obtained

Publisher bot walls blocked the full text of Tsirogiannis, Sandel and Cheliotis
(2012), Tsirogiannis, Sandel and Kalvisa (2014), Tsirogiannis and Sandel (2016,
RECOMB; the non-uniform null) and Tsirogiannis and Sandel (2016, Ecography; the
package paper). The 2014 algorithm is known here only through the description in
the open-access PLoS ONE paper and through the source code. Nothing above
depends on them, since the formulas are verified against the package itself, but
the 2014 paper would be worth reading before optimizing.

## Sources

* Sandel, B. (2018). Richness-dependence of phylogenetic diversity indices.
  Ecography 41: 837-844. doi:10.1111/ecog.02967
* Nipperess, D. A. and Matsen, F. A. (2013). The mean and variance of
  phylogenetic diversity under rarefaction. Methods in Ecology and Evolution 4:
  566-572. doi:10.1111/2041-210X.12042; arXiv:1208.6552
* Tsirogiannis, C. and Sandel, B. (2016). Fast computations for measures of
  phylogenetic beta diversity. PLoS ONE 11: e0151167.
  doi:10.1371/journal.pone.0151167
* Tsirogiannis, C. and Sandel, B. (2014). Computing the skewness of the
  phylogenetic mean pairwise distance in linear time. Algorithms for Molecular
  Biology 9: 15. doi:10.1186/1748-7188-9-15
* Tsirogiannis, C., Sandel, B. and Kalvisa, A. (2014). New algorithms for
  computing phylogenetic biodiversity. WABI 2014, LNCS 8701: 187-203.
  doi:10.1007/978-3-662-44753-6_15
* Tsirogiannis, C., Sandel, B. and Cheliotis, D. (2012). Efficient computation
  of popular phylogenetic tree measures. WABI 2012, LNCS 7534: 30-43.
  doi:10.1007/978-3-642-33122-0_3
* Tsirogiannis, C. and Sandel, B. (2016). PhyloMeasures: a package for computing
  phylogenetic biodiversity measures and their statistical moments. Ecography
  39: 709-714. doi:10.1111/ecog.01814
* PhyloMeasures 2.1 source, CRAN archive.
* picante 1.8.4, documentation for `expected.pd`.
