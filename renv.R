#!/usr/bin/env Rscript
# install_packages.sh
# This script sets up a reproducible R environment with renv
# and installs the required packages.

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialize renv (creates a new environment and renv.lock file)
renv::init(bare = TRUE)

# Define the CRAN packages to install
direct_packages <- c(
  "glue", "ggplot2", "ggpubr", "jsonlite", "optparse", "pbmcapply", "readr"
)

indirect_packages <- c(
  "askpass",
  "base64enc",
  "brew",
  "brio",
  "bslib",
  "cachem",
  "callr",
  "cli",
  "clipr",
  "commonmark",
  "cpp11",
  "crayon",
  "credentials",
  "curl",
  "desc",
  "devtools",
  "diffobj",
  "digest",
  "downlit",
  "ellipsis",
  "evaluate",
  "fansi",
  "fastmap",
  "fontawesome",
  "fs",
  "gert",
  "gh",
  "gitcreds",
  "glue",
  "highr",
  "htmltools",
  "htmlwidgets",
  "httpuv",
  "httr2",
  "ini",
  "jquerylib",
  "jsonlite",
  "knitr",
  "later",
  "lifecycle",
  "magrittr",
  "memoise",
  "mime",
  "miniUI",
  "openssl",
  "pillar",
  "pkgbuild",
  "pkgconfig",
  "pkgdown",
  "pkgload",
  "praise",
  "prettyunits",
  "processx",
  "profvis",
  "promises",
  "ps",
  "purrr",
  "R6",
  "ragg",
  "rappdirs",
  "rcmdcheck",
  "Rcpp",
  "remotes",
  "rlang",
  "rmarkdown",
  "roxygen2",
  "rprojroot",
  "rstudioapi",
  "rversions",
  "sass",
  "sessioninfo",
  "shiny",
  "sourcetools",
  "stringi",
  "stringr",
  "sys",
  "systemfonts",
  "testthat",
  "textshaping",
  "tibble",
  "tinytex",
  "urlchecker",
  "usethis",
  "utf8",
  "vctrs",
  "waldo",
  "whisker",
  "withr",
  "xfun",
  "xml2",
  "xopen",
  "xtable",
  "yaml",
  "zip"
)

cran_packages <- c(
  direct_packages,
  indirect_packages
)

# Install the CRAN packages
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# Ensure devtools is installed for GitHub installations
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Install the GitHub package from mrc-ide
devtools::install_github("YuLab-SMU/ggtree@c17773c973d6c4036ee3af40a3957fb74d8ee9ff")
devtools::install_github('mmore500/mlscluster@ddddeeb0e26c26c6f5371de16b06b8e4c88c26a4')

# Snapshot the current state of the library to record package versions in renv.lock
renv::snapshot(packages = cran_packages)
