# please install those packages if you don't have them
libs_load <- c("mlscluster", "glue", "ggplot2", "ggpubr")
invisible(lapply(libs_load, library, character.only = TRUE))

NCPU <- 7
PERIOD_INTEREST <- 2  # somehow needed in stats_multiple_thresholds, wtf
THRESHOLD_INTEREST <- 2
options(scipen = 999)

thr <- c(0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
quantl <- c(
  1 / 400,
  1 / 200,
  1 / 400,
  1 / 100,
  2 / 100,
  3 / 100,
  4 / 100,
  5 / 100,
  1 / 10,
  1 / 4
)
pal_lineages <- c(
  B.1.1.7 = "#fd7f6f",
  `AY.4.*` = "#7eb0d5",
  `AY.* (non-AY.4.*)` = "#b2e061",
  `B.1.177.*` = "#bd7ebe",
  `BA.1.*` = "#ffb55a",
  `BA.2.*` = "#ffcc00",
  Other = "#fdcce5"
)
thr_pref <- "threshold_quantile"
path_thresholds <- c(
  glue::glue("{thr_pref}_0.25"),
  glue::glue("{thr_pref}_0.5"),
  glue::glue("{thr_pref}_0.75"),
  glue::glue("{thr_pref}_1"),
  glue::glue("{thr_pref}_2"),
  glue::glue("{thr_pref}_3"),
  glue::glue("{thr_pref}_4"),
  glue::glue("{thr_pref}_5"),
  glue::glue("{thr_pref}_10"),
  glue::glue("{thr_pref}_25")
)

# Load metadata and tree sc2_md_curated <-
# readRDS('rds/sc2_md_curated_WITH_Xs_Ns.rds')
# sc2_tre_curated <- readRDS('rds/sc2_tre_curated.rds')
sc2_md_curated <- readRDS(url(
  "https://osf.io/f953r/download",
  "rb"
))
sc2_tre_curated <- readRDS(url(
  "https://osf.io/24r3e/download",
  "rb"
))

# options(max.print=999999)

## 1. Run mlsclust ## Period 1: testing (without artifact
## removal: ~10 minutes: 6 for mlsclust and 4 for
## run_diff_thresholds)
start <- Sys.time()
res_p1 <- mlsclust(
  sc2_tre_curated,
  sc2_md_curated,
  min_descendants = 10,
  max_descendants = 20000,
  min_cluster_age_yrs = 1 / 12,
  min_date = as.Date("2019-12-30"),
  max_date = as.Date("2020-12-31"),
  branch_length_unit = "days",
  rm_seq_artifacts = TRUE,
  defining_mut_threshold = 0.75,
  root_on_tip = "Wuhan/WH04/2020",
  root_on_tip_sample_time = 2019.995,
  detailed_output = FALSE,
  ncpu = NCPU
) # detailed_output=TRUE
end <- Sys.time()
total_time_p1 <- as.numeric(end - start, units = "mins")
message(paste("Total time elapsed:", total_time_p1, "mins"))
system(glue("mkdir -p rds/"))
saveRDS(res_p1, "rds/res_p1.rds")

# res_p1 <- readRDS('rds/res_p1.rds')

start <- Sys.time()
for (i in 1:length(thr)) {
  run_diff_thresholds(
    stats_df_unfilt = res_p1[[1]],
    tgt_nodes = res_p1[[2]],
    homoplasies = res_p1[[3]],
    output_dir = glue(
      "results/period1/threshold_quantile_{thr[[i]]}/"
    ),
    quantile_choice = quantl[i],
    quantile_threshold_ratio_sizes = thr[i],
    quantile_threshold_ratio_persist_time = thr[i],
    quantile_threshold_logit_growth = thr[i],
    threshold_keep_lower = TRUE,
    compute_tree_annots = FALSE,
    plot_global_mut_freq = FALSE,
    detailed_output = FALSE
  ) # , desc_sisters=list(res_p1[[4]], res_p1[[5]]), amd=res_p1[[6]]
}
end <- Sys.time()
total_time_p1_t <- as.numeric(end - start, units = "mins")
message(paste(
  "Total time elapsed:",
  total_time_p1_t,
  "mins"
))

## 2. Run statistical tests for multiple thresholds ##
# Expected paths containing results for each quantile threshold
table_names_periods <- c(2)
table_names_thresholds <- c(
  0.25,
  0.5,
  0.75,
  1,
  2,
  3,
  4,
  5,
  10,
  25
)
table_combs <- do.call(
  paste,
  c(
    expand.grid(
      table_names_periods,
      table_names_thresholds
    ),
    sep = "_"
  )
)
out_folder <- "period1"

system(glue("mkdir -p stat_results/plots_paper/"))
# Override pal_lineages from package with new major_lineage labels
major_lineages <- c(  # somehow needed in stats_multiple_thresholds, wtf
  "B.1.177.*",
  "B.1.1.7",
  "AY.4.*",
  "AY.* (non-AY.4.*)",
  "Other"
) # 'Beta_B.1.351','Gamma_P.1'
rmult_p2 <- stats_multiple_thresholds(
  "results/period1",
  pal_lineages,
  "period1",
  "ED_FileS1.txt"
)
plot_tbc_stats("results/period1", out_folder)
stacked_nsites_genomic_region_mult_thresholds(
  "stat_results/period1/genomewide_plot_non-syn/"
)

## 4. Run for the selected threshold (2% since FDR ~10% and epsilon ~2.5%)
# Expected paths containing results for each quantile threshold
out_folder <- "period1_thr2"

r2 <- stats_selected_threshold(
  "results/period1",
  thr_index = 5,
  pal_lineages,
  out_folder
) # NOTE: thr_index = 5 -> 2%
# Extended data fig. S5: mlscluster, intersection, and
# hyphy (also stratified by lineage) identified sites for
# thr=2% and period before Omicron
r2[[4]]
ggsave(
  glue("stat_results/plots_paper/ED_FigS5.png"),
  units = "px",
  width = 2200,
  height = 2500,
  dpi = 300,
  bg = "white"
)
