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
cran_packages <- c("glue", "ggplot2", "ggpubr", "pbmcapply")

# Install the CRAN packages
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# Ensure devtools is installed for GitHub installations
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Install the GitHub package from mrc-ide
devtools::install_github("YuLab-SMU/ggtree@c17773c973d6c4036ee3af40a3957fb74d8ee9ff")
devtools::install_github('mrc-ide/mlscluster@e9f7939ab69a61bc310a968252d96177c59921ba')

# Snapshot the current state of the library to record package versions in renv.lock
renv::snapshot()
