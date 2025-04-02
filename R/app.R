#!/usr/bin/env Rscript

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000) # for gzcon compat, wtf
timeout <- 600
max_tries <- 10
options(timeout = timeout)

library("optparse")

# Define command line options for mlsclust arguments and input data sources.
option_list = list(
  make_option(
    c("--metadata_json"),
    type = "character",
    default = "gzcon(url('https://osf.io/xvh7d/download', 'rb'))",
    help = "File path or url for metadata (JSON file) [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--phylogeny_json"),
    type = "character",
    default = "gzcon(url('https://osf.io/85s3q/download', 'rb'))",
    help = "File path or url for phylogeny tree (JSON file) [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--min_descendants"),
    type = "character",
    default = "10",
    help = "Minimum number of descendants [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--max_descendants"),
    type = "character",
    default = "20000",
    help = "Maximum number of descendants [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--min_cluster_age_yrs"),
    type = "character",
    default = "1/12",
    help = "Minimum cluster age in years [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--min_date"),
    type = "character",
    default = "as.Date('2019-12-30')",
    help = "Minimum date [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--max_date"),
    type = "character",
    default = "as.Date('2020-12-31')",
    help = "Maximum date [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--branch_length_unit"),
    type = "character",
    default = "'days'",
    help = "Branch length unit [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--rm_seq_artifacts"),
    type = "character",
    default = "TRUE",
    help = "Remove sequence artifacts? [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--defining_mut_threshold"),
    type = "character",
    default = "0.75",
    help = "Defining mutation threshold [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--root_on_tip"),
    type = "character",
    default = "'Wuhan/WH04/2020'",
    help = "Tip to root on [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--root_on_tip_sample_time"),
    type = "character",
    default = "2019.995",
    help = "Sample time for the root tip [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--detailed_output"),
    type = "character",
    default = "FALSE",
    help = "Detailed output? [default: %default]",
    metavar = "character"
  ),
  make_option(
    c("--ncpu"),
    type = "character",
    default = "NA",
    help = "Number of CPU cores [default: detected automatically]",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Evaluate each argument so that types are correctly set.
metadata_json <- eval(parse(text = opt$metadata_json))
phylogeny_json <- eval(parse(text = opt$phylogeny_json))
min_descendants <- eval(parse(text = opt$min_descendants))
max_descendants <- eval(parse(text = opt$max_descendants))
min_cluster_age_yrs <- eval(parse(text = opt$min_cluster_age_yrs))
min_date <- eval(parse(text = opt$min_date))
max_date <- eval(parse(text = opt$max_date))
branch_length_unit <- eval(parse(text = opt$branch_length_unit))
rm_seq_artifacts <- eval(parse(text = opt$rm_seq_artifacts))
defining_mut_threshold <- eval(parse(text = opt$defining_mut_threshold))
root_on_tip <- eval(parse(text = opt$root_on_tip))
root_on_tip_sample_time <- eval(parse(text = opt$root_on_tip_sample_time))
detailed_output <- eval(parse(text = opt$detailed_output))

if (opt$ncpu == "NA") {
  library(parallel)
  NCPU <- detectCores()
  if (is.na(NCPU)) {
    NCPU <- 7 # Fallback if detection fails
  }
  ncpu <- NCPU
} else {
  ncpu <- eval(parse(text = opt$ncpu))
}

print(
  list(
    metadata_json = metadata_json,
    phylogeny_json = phylogeny_json,
    min_descendants = min_descendants,
    max_descendants = max_descendants,
    min_cluster_age_yrs = min_cluster_age_yrs,
    min_date = min_date,
    max_date = max_date,
    branch_length_unit = branch_length_unit,
    rm_seq_artifacts = rm_seq_artifacts,
    defining_mut_threshold = defining_mut_threshold,
    root_on_tip = root_on_tip,
    root_on_tip_sample_time = root_on_tip_sample_time,
    detailed_output = detailed_output
  )
)

###############################################################################
# Load required libraries
###############################################################################
libs_load <- c(
  "mlscluster",
  "glue",
  "ggplot2",
  "ggpubr",
  "jsonlite",
  "readr"
)
invisible(lapply(libs_load, library, character.only = TRUE))

options(scipen = 999)

###############################################################################
message(">>> Run mlsclust")
###############################################################################

# Load metadata and tree using the provided URLs (or file paths)
for (i in 1:max_tries) {
  tryCatch(
    {
      metadata_df <- unserializeJSON(read_lines(eval(parse(
        text = opt$metadata_json
      ))))
      break
    },
    error = function(e) {
      message("Error reading metadata JSON file: ", e)
      if (i == max_tries) {
        stop("Max tries reached. Exiting.")
      }
    }
  )
}
system(glue("mkdir -p rds/"))
writeLines(serializeJSON(metadata_df), "rds/metadata_df.json")
saveRDS(metadata_df, "rds/metadata_df.rds")

for (i in 1:max_tries) {
  tryCatch(
    {
      phylogeny_tree <- unserializeJSON(read_lines(eval(parse(
        text = opt$phylogeny_json
      ))))
      break
    },
    error = function(e) {
      message("Error reading phylogeny JSON file: ", e)
      if (i == max_tries) {
        stop("Max tries reached. Exiting.")
      }
    }
  )
}
system(glue("mkdir -p rds/"))
writeLines(serializeJSON(phylogeny_tree), "rds/phylogeny_tree.json")
saveRDS(phylogeny_tree, "rds/phylogeny_tree.rds")

# ----------------------------------------------------------------------------

start <- Sys.time()
mlsclust_result <- mlsclust(
  phylogeny_tree,
  metadata_df,
  min_descendants = min_descendants,
  max_descendants = max_descendants,
  min_cluster_age_yrs = min_cluster_age_yrs,
  min_date = min_date,
  max_date = max_date,
  branch_length_unit = branch_length_unit,
  rm_seq_artifacts = rm_seq_artifacts,
  defining_mut_threshold = defining_mut_threshold,
  root_on_tip = root_on_tip,
  root_on_tip_sample_time = root_on_tip_sample_time,
  detailed_output = detailed_output,
  ncpu = ncpu
)
# mlsclust_result <- readRDS('rds/mlsclust_result.rds')
message(paste(
  "Total time elapsed:",
  as.numeric(Sys.time() - start, units = "mins"),
  "mins"
))
system(glue("mkdir -p rds/"))
saveRDS(mlsclust_result, "rds/mlsclust_result.rds")
writeLines(serializeJSON(mlsclust_result), "rds/mlsclust_result.json")

# test: pass through JSON
mlsclust_result <- unserializeJSON(serializeJSON(mlsclust_result))

# unpack mlsclust_result
clustering_statistics_df <- mlsclust_result[[1]]
target_nodes_vector <- mlsclust_result[[2]]
homoplasy_freq_df <- mlsclust_result[[3]]

writeLines(
  serializeJSON(clustering_statistics_df),
  "rds/clustering_statistics_df.json"
)
writeLines(
  serializeJSON(target_nodes_vector),
  "rds/target_nodes_vector.json"
)
writeLines(
  serializeJSON(homoplasy_freq_df),
  "rds/homoplasy_freq_df.json"
)

###############################################################################
message(">>> Thresholds")
###############################################################################
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

# ----------------------------------------------------------------------------

start <- Sys.time()
for (i in 1:length(thr)) {
  run_diff_thresholds(
    stats_df_unfilt = clustering_statistics_df,
    tgt_nodes = target_nodes_vector,
    homoplasies = homoplasy_freq_df,
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
  )
}
message(paste(
  "Total time elapsed:",
  as.numeric(Sys.time() - start, units = "mins"),
  "mins"
))

###############################################################################
message(">>> Run statistical tests for multiple thresholds")
###############################################################################

# somehow needed in stats_multiple_thresholds, wtf
major_lineages <- c(
  "B.1.177.*",
  "B.1.1.7",
  "AY.4.*",
  "AY.* (non-AY.4.*)",
  "Other"
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

# Expected paths containing results for each quantile threshold
table_names_periods <- c(2)
table_names_thresholds <- thr
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

PERIOD_INTEREST <- 2 # somehow needed in stats_multiple_thresholds, wtf

# ----------------------------------------------------------------------------

system(glue("mkdir -p stat_results/plots_paper/"))
_ <- stats_multiple_thresholds(
  "results/period1",
  pal_lineages, # Override pal_lineages from package with new labels
  "period1",
  "ED_FileS1.txt"
)
plot_tbc_stats("results/period1", "period1")
stacked_nsites_genomic_region_mult_thresholds(
  "stat_results/period1/genomewide_plot_non-syn/"
)


###############################################################################
message(">>> Selected threshold (2% since FDR ~10% & epsilon ~2.5%)")
###############################################################################

THRESHOLD_INTEREST <- 2

# ----------------------------------------------------------------------------

_ <- stats_selected_threshold(
  "results/period1",
  thr_index = 5,
  pal_lineages,
  "period1_thr2"
) # NOTE: thr_index = 5 -> 2%

# Extended data fig. S5: mlscluster, intersection, and hyphy (also stratified
# by lineage) identified sites for thr=2%
ggsave(
  glue("stat_results/plots_paper/ED_FigS5.png"),
  units = "px",
  width = 2200,
  height = 2500,
  dpi = 300,
  bg = "white"
)
