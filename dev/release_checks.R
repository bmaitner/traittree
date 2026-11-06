# Pre-release checklist for traittree. Source the sections one at a time from
# the package root. This folder is in .Rbuildignore and never ships.
#
# The GitHub Actions workflow (.github/workflows/R-CMD-check.yaml) already runs
# R CMD check on Linux/macOS/Windows for every push, so the steps here are the
# ones CRAN asks for beyond that: the --as-cran incoming checks, URL and
# spelling checks, tarball size, and the remote CRAN-like builders.

# 0. The archived dependency ------------------------------------------------
# PhyloMeasures provides the phylogenetic diversity queries. It was archived
# from CRAN in 2022 and is a *suggested* dependency, so the package installs
# and checks without it and calculate_pd_metric_rasters() errors with install
# instructions. Local checks should have it present, so that the diversity
# code is actually exercised:
# remotes::install_url("https://cran.r-project.org/src/contrib/Archive/PhyloMeasures/PhyloMeasures_2.1.tar.gz")
#
# This is the one thing standing between traittree and CRAN. CRAN requires
# suggested packages to be installable from a declared repository; an archived
# one is not, so a submission would need either Additional_repositories, or
# PhyloMeasures reinstated on CRAN, or those two pd.query() calls moved to a
# maintained package. Settle that before spending effort on the steps below.

# 1. Bookkeeping ------------------------------------------------------------
# Bump Version in DESCRIPTION and add a NEWS.md entry. There is no
# cran-comments.md yet; usethis::use_cran_comments() creates one, and it is
# where the PhyloMeasures situation has to be explained to CRAN.
# usethis::use_version("minor")        # or "patch" / "major"
# usethis::use_release_issue()         # tickable checklist as a GitHub issue

# 2. Regenerate docs and data ------------------------------------------------
devtools::document()
source("data-raw/generate_example_data.R")   # example_tree / traits / occurrences
# There is no pkgdown site; usethis::use_pkgdown_github_pages() sets one up.

# 3. Local checks -----------------------------------------------------------
devtools::test()
chk <- devtools::check(cran = TRUE, incoming = TRUE, remote = TRUE)
stopifnot(length(chk$errors) == 0, length(chk$warnings) == 0)

# The suggested dependency has to be genuinely optional, not just documented as
# optional. This is the check that catches a pd.query() call escaping its
# guard, or an example or vignette chunk that forgets to skip. It hides every
# Suggests, so testthat and rmarkdown go too: expect the tests and the vignette
# to be skipped rather than to pass, and read the failures rather than the
# status. For a version that hides PhyloMeasures alone, build a library of
# junctions to every other installed package and point R_LIBS_USER at it.
# Sys.setenv("_R_CHECK_DEPENDS_ONLY_" = "true"); devtools::check()

# URLs in Rd files, README, and DESCRIPTION. The CRAN archive link for
# PhyloMeasures is expected to resolve; if it ever stops, the install
# instructions in the README and in require_phylomeasures() are both wrong.
urlchecker::url_check()

# Spelling. There is no inst/WORDLIST yet, so the first run reports every
# domain term (phylogeny, ultrametric, Brownian, phylopars, ...). Create the
# list, then keep it honest.
spelling::spell_check_package()
# spelling::update_wordlist()          # accept the current list

# Tarball size. Nothing large ships: the example data are a few KB and there
# are no figures, so this is slack unless something unexpected creeps in.
tb <- devtools::build(path = tempdir(), quiet = TRUE)
cat(basename(tb), round(file.size(tb) / 1024), "KB\n")
stopifnot(file.size(tb) < 4 * 1024^2)

# The vignette does the real end-to-end run. Worth reading the rendered output
# rather than only checking that it builds, since the numbers it quotes in the
# prose (correlations, uncertainty) come from the data and can drift.
# devtools::build_vignettes(); browseURL("doc/traittree.html")

# 4. Remote builders (results arrive by email / GitHub) ---------------------
devtools::check_win_devel()            # CRAN's Windows R-devel builder

# R-hub v2 runs in this repo's GitHub Actions (.github/workflows/rhub.yaml,
# created once with rhub::rhub_setup()). The workflow file must be on the
# default branch, so commit and push before checking. Results appear under
# https://github.com/bmaitner/traittree/actions.
#
# "nosuggests" is the one to watch here: it is the hosted version of the
# depends-only check above, and the reason this package can afford to keep an
# archived dependency at all.
rhub::rhub_doctor()
rhub::rhub_check(platforms = c("linux", "macos", "windows",
                               "ubuntu-next", "nosuggests", "donttest"))
# rhub::rhub_platforms()               # the full list of available builders

# 5. Submit -----------------------------------------------------------------
# devtools::release()                  # interactive checklist, then uploads
# After acceptance: usethis::use_github_release(); usethis::use_dev_version()
