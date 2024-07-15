# mlscluster-experiments

[![DOI](https://zenodo.org/badge/685186608.svg)](https://zenodo.org/doi/10.5281/zenodo.10520059)

This repository contains the code needed to reproduce the results in the *mlscluster* paper.

Fork this repo and install the `mlscluster` package using `devtools::install_github("mrc-ide/mlscluster")` to have everything set up.

Also download "Extended Data File S3" (zip file) from [this Zenodo link](https://zenodo.org/doi/10.5281/zenodo.10276239) and place it within a `rds` folder within your recently forked repo. `sc2_md_curated_WITH_Xs_Ns.rds` contains necessary metadata (columns: sequence_name, sample_date_lineage, major_lineage, mutations). `sc2_tre_curated.rds` contains the time-scaled tree estimated for these sequences. `res_p2.rds` and `res_p3.rds` contain outputs of `mlsclust` for period excluding Omicron (June 2020 to mid-November 2021) and including it until Pillar2 termination (June 2020 to April 2022), respectively. See **NOTE** below for more details.

Figures 3 to 5 and Extended Data outputs can be replicated following the `R/running.R` script.

Code for comparisons against deep mutational scanning experiments and comprehensive between-host fitness estimates from [this paper](https://academic.oup.com/ve/article/9/2/vead055/7265011) can be found in `R/dms_fitness_comparisons.R`.

Although the script `R/spatiotemporal_distr_samples.R` has the code to reproduce Fig. 2, it uses non-public adm2-level geographic data, which is only accessible for COG-UK registered members.

**NOTE**: 

* `R/annotate_sites_missing_calls.R` was previously used to perform the alignment-aware artifact removal method, so that `sc2_md_curated_WITH_Xs_Ns.rds` already includes "X" and "N" missing sites within its `mutations` column (see paper "Methods: Statistical analysis for identifying genomic regions enriched for TFPs") and is ready for proper analysis
* We successfully run the `mlsclust` function for ~66 thousand sequences (June 2020 to December 2020) using a laptop with core i7, 32 GB of RAM and 8 CPU cores in ~15 minutes. <br>
Since the dataset with >600 thousand sequences (June 2020 to mid-November 2021, before Omicron) and the dataset comprising >1.2 million sequences (June 2020 to April 2022, including Omicron) took ~15 and ~40 hours respectively in this setting, we distributed the `res_p2.rds` and `res_p3.rds` outputs to facilitate reproducibility. To replicate and obtain these outputs on your own, we recommend using a server with higher computational power (e.g. 256 GB of RAM and 16 CPUs).