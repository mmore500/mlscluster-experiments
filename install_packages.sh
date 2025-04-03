#!/bin/bash

# install_packages.sh
# Remove renv and install all R packages globally

set -euo pipefail  # Exit on error

echo "Removing renv files..."
rm -rf renv/ renv.lock .Rprofile

echo "Installing R packages globally..."

Rscript -e '
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
'

Rscript -e '
# Define CRAN packages
direct_packages <- c(
  "glue", "ggplot2", "ggpubr", "jsonlite", "optparse", "pbmcapply", "readr"
)

indirect_packages <- c(
  "askpass", "base64enc", "brew", "brio", "bslib", "cachem", "callr", "cli",
  "clipr", "commonmark", "cpp11", "crayon", "credentials", "curl", "desc",
  "diffobj", "digest", "downlit", "ellipsis", "evaluate", "fansi",
  "fastmap", "fontawesome", "fs", "gert", "gh", "gitcreds", "glue", "highr",
  "htmltools", "htmlwidgets", "httpuv", "httr2", "ini", "jquerylib", "jsonlite",
  "knitr", "later", "lifecycle", "magrittr", "memoise", "mime", "miniUI",
  "openssl", "pillar", "pkgbuild", "pkgconfig", "pkgdown", "pkgload", "praise",
  "prettyunits", "processx", "profvis", "promises", "ps", "purrr", "R6", "ragg",
  "rappdirs", "rcmdcheck", "Rcpp", "remotes", "rlang", "rmarkdown", "roxygen2",
  "rprojroot", "rstudioapi", "rversions", "sass", "sessioninfo", "shiny",
  "sourcetools", "stringi", "stringr", "sys", "systemfonts", "testthat",
  "textshaping", "tibble", "tinytex", "urlchecker", "usethis", "utf8", "vctrs",
  "waldo", "whisker", "withr", "xfun", "xml2", "xopen", "xtable", "yaml", "zip"
)

cran_packages <- c(direct_packages, indirect_packages)

# Install CRAN packages
install.packages(cran_packages, repos = "https://cloud.r-project.org")
'

Rscript -e '
# Install GitHub packages
devtools::install_github("YuLab-SMU/ggtree@c17773c973d6c4036ee3af40a3957fb74d8ee9ff")
devtools::install_github("mmore500/mlscluster@85a39581ed84726c29afb8b5ae74b5f524b998df")
'
