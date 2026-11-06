# Research check: Sandel (2018) removes richness-dependence by rarefying every
# community to a common richness k and averaging the index over replicates (his
# "RR", 1000 replicates). For fixed k the index is an affine function of PD, so
# the average over rarefactions should equal
#     ( E_k[PD of the community's own subtree] - mu_k(T) ) / sigma_k(T)
# with the numerator's expectation given exactly by the same formula applied to
# the community subtree. Check that against brute-force replication.
suppressMessages(library(ape))
load("data/example_tree.rda"); tree <- example_tree
s <- length(tree$tip.label)

expected_pd <- function(tr, r) {            # Nipperess & Matsen eq. 5 (unrooted)
  n <- length(tr$tip.label)
  nd <- ape::node.depth(tr)[tr$edge[, 2]]    # tips below each edge
  h <- function(x) ifelse(x < r, 0, exp(lchoose(x, r) - lchoose(n, r)))
  sum(tr$edge.length * (1 - h(nd) - h(n - nd)))
}
pd <- function(tr, spp) if (length(spp) < 2) 0 else sum(ape::keep.tip(tr, spp)$edge.length)

set.seed(11); k <- 10
mom <- PhyloMeasures::pd.moments(tree, sample.sizes = k)
cat(sprintf("%8s %12s %12s %12s\n", "richness", "closed form", "RR (4000)", "difference"))
for (S in c(15, 30, 45)) {
  community <- sample(tree$tip.label, S)
  sub <- ape::keep.tip(tree, community)
  exact <- (expected_pd(sub, k) - mom[1, 1]) / mom[1, 2]
  reps <- replicate(4000, (pd(tree, sample(community, k)) - mom[1, 1]) / mom[1, 2])
  cat(sprintf("%8d %12.5f %12.5f %12.5f   (MC s.e. %.4f)\n", S, exact, mean(reps), exact - mean(reps), sd(reps)/sqrt(4000)))
}
