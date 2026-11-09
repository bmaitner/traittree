# Pre-release checklist for traittree. Source the sections one at a time from
# the package root. This folder is in .Rbuildignore and never ships.
#
# The GitHub Actions workflow (.github/workflows/R-CMD-check.yaml) already runs
# R CMD check on Linux/macOS/Windows for every push, so the steps here are the
# ones CRAN asks for beyond that: the --as-cran incoming checks, URL and
# spelling checks, tarball size, and the remote CRAN-like builders.

# 1. Bookkeeping ------------------------------------------------------------
# Bump Version in DESCRIPTION and add a NEWS.md entry, then bring
# cran-comments.md up to date with the results of the checks below.
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

# URLs in Rd files, README, and DESCRIPTION.
urlchecker::url_check()

# Spelling. There is no inst/WORDLIST yet, so the first run reports every
# domain term (phylogeny, ultrametric, Brownian, phylopars, ...). Create the
# list, then keep it honest.
spelling::spell_check_package()
# spelling::update_wordlist()          # accept the current list

# The phylogenetic diversity code is tested against values frozen from
# PhyloMeasures 2.1 (tests/testthat/fixtures/pd_reference.rds). Those never need
# regenerating unless the fixtures themselves change;
# dev/pdi_research/make_pd_fixtures.R says how, and dev/pdi_research_notes.md
# records where the formulas come from.

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
# These cover what the GitHub Actions matrix cannot: CRAN's own machines, and
# builds whose arithmetic differs. Everything here checks what is on GitHub or
# in the working tree, so commit and push first.

# CRAN's Windows machines, on all three R versions, and its macOS builder.
devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()
devtools::check_mac_release()
devtools::check_mac_devel()

# win-builder also runs the spell check of DESCRIPTION, which a Windows machine
# without aspell skips. The surnames and "Phylogenies" it would flag are in
# .aspell/traittree.rds (see .aspell/defaults.R for how to change the list).

# R-hub v2 runs in this repo's GitHub Actions (.github/workflows/rhub.yaml,
# created once with rhub::rhub_setup()). The workflow file must be on the
# default branch. Results appear under
# https://github.com/bmaitner/traittree/actions.
#
# Beyond the three operating systems on R-devel, the platforms are chosen for
# numerics, because the tests compare against frozen values at tight tolerances
# and run an optimizer (Rphylopars): atlas and mkl swap the linear algebra
# library, nold has no long doubles. CRAN checks on all three. The memory
# checkers (clang-asan, valgrind, rchk) are left out: there is no compiled code.
#
# The Fedora containers (atlas, mkl, nold, nosuggests) can fail before the
# package is built, when a prebuilt dependency does not match the container's
# R: in November 2026 nosuggests died loading ape ("libRlapack.so: cannot open
# shared object file"). That says nothing about the package; a failure only
# counts if the log gets as far as "checking for file 'traittree/DESCRIPTION'".
rhub::rhub_doctor()
rhub::rhub_check(platforms = c("linux", "macos", "windows",
                               "ubuntu-next", "donttest",
                               "atlas", "mkl", "nold"))
# rhub::rhub_platforms()               # the full list of available builders

# 5. Submit -----------------------------------------------------------------
# devtools::release()                  # interactive checklist, then uploads
# After acceptance: usethis::use_github_release(); usethis::use_dev_version()
