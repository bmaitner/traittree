# Research check, not package code: evaluate the published formulas directly
# (Nipperess & Matsen 2013, eq. 5 and Appendix) by brute force over edge pairs,
# and compare with PhyloMeasures::pd.moments and with Monte Carlo.
suppressMessages(library(ape))
load("data/example_tree.rda"); tree <- example_tree
s <- length(tree$tip.label); w <- tree$edge.length; ne <- nrow(tree$edge)

# D[e, t] = 1 if tip t descends from edge e
D <- matrix(0, ne, s)
for (t in 1:s) { path <- ape::nodepath(tree, from = s + 1, to = t); D[match(path[-1], tree$edge[, 2]), t] <- 1 }
se <- rowSums(D); I <- D %*% t(D)                      # clade sizes, pairwise overlaps

moments <- function(r) {
  h <- function(x) ifelse(x < r, 0, exp(lchoose(pmax(x, 0), r) - lchoose(s, r)))   # P(all r draws fall in a set of x tips)
  p <- 1 - h(se) - h(s - se)                                                         # P(edge in spanning subtree)
  mu <- sum(w * p)
  Se <- matrix(se, ne, ne); Sf <- t(Se)
  nested_e_in_f <- (I == Se) & (Se < Sf | (Se == Sf & row(I) != col(I)))             # e strictly below f (or twin)
  both <- 1 - h(s - Se) - h(s - Sf) + h(s - Se - Sf)                                 # disjoint clades
  lo <- pmin(Se, Sf); hi <- pmax(Se, Sf)
  nested <- (I == lo) & (I > 0)
  both[nested] <- (1 - h(hi) - h(s - lo) + h(hi - lo))[nested]                       # one clade inside the other
  diag(both) <- p
  c(mean = mu, sd = sqrt(max(sum(outer(w, w) * both) - mu^2, 0)))
}

rs <- c(2, 3, 8, 20, 45, 59, 60)
mine <- t(sapply(rs, moments))
pm <- PhyloMeasures::pd.moments(tree, sample.sizes = rs)
cat(sprintf("%3s  %10s %10s  %10s %10s\n", "r", "mean(mine)", "mean(PM)", "sd(mine)", "sd(PM)"))
for (i in seq_along(rs)) cat(sprintf("%3d  %10.6f %10.6f  %10.6f %10.6f\n", rs[i], mine[i,1], pm[i,1], mine[i,2], pm[i,2]))
cat("\nmax |mean diff|:", signif(max(abs(mine[,1] - pm[,1])), 3), "  max |sd diff|:", signif(max(abs(mine[,2] - pm[,2])), 3), "\n")

# independent Monte Carlo confirmation at r = 8
set.seed(1); r <- 8
mc <- replicate(20000, { tips <- sample(s, r); sum(w[rowSums(D[, tips, drop = FALSE]) > 0 & rowSums(D[, tips, drop = FALSE]) < r]) })
cat(sprintf("Monte Carlo r=8 (20000 draws): mean %.4f (exact %.4f)   sd %.4f (exact %.4f)\n", mean(mc), mine[3,1], sd(mc), mine[3,2]))
