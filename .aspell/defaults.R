# Read by R's spell check of DESCRIPTION, including the CRAN incoming check
# (utils:::aspell_package_description). "en_stats" is the dictionary that check
# uses by default; "traittree" is .aspell/traittree.rds, the words below.
# To change the list, edit the words and run this from the package root:
#   saveRDS(c("Cadotte", "Matsen", "Nipperess", "Phylogenies", "Sandel"),
#           ".aspell/traittree.rds", version = 2)
description <- list(dictionaries = c("en_stats", "traittree"))
