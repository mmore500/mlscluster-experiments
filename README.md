# mlscluster-experiments

This repository contains the code and data needed to reproduce the results in the *mlscluster* paper.

Fork this repo and install the mlscluster package using `devtools::install_github("vinibfranc/mlscluster")` to have everything set up.

The `R/spatiotemporal_distr_samples.R` script reproduces Figure 2.

Figures 3 to 5 and supplementary information can be replicated following the `R/running.R` script.

**Note**: We successfully run the `mlsclust` function for ~66 thousand sequences (June 2020 to December 2020) using a laptop with core i7, 32 GB of RAM and 8 CPU cores (4 used) in ~6 minutes. <br>
Since the dataset with >600 thousand sequences (June 2020 to mid-November 2021, before Omicron) and the dataset comprising >1.2 million sequences (June 2020 to April 2022, including Omicron) took ~10 and ~30 hours respectively in this setting, we distributed the `res_p2` and `res_p3` outputs within the *mlscluster* package to facilitate analysis. To replicate and obtain these outputs on your own, we recommend using a server with higher computational power (e.g. 96 GB of RAM and 16 CPUs).