libs_load <- c("mlscluster","glue", "ggplot2","ggpubr") # please install those packages if you don't have them
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 7
options(scipen=999)

thr <- c(0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
quantl <- c(1/400, 1/200, 1/400, 1/100, 2/100, 3/100, 4/100, 5/100, 1/10, 1/4)
pal_lineages <- c("B.1.1.7"="#fd7f6f", "AY.4.*"="#7eb0d5", "AY.* (non-AY.4.*)"="#b2e061", "B.1.177.*"="#bd7ebe", "BA.1.*"="#ffb55a", "BA.2.*"="#ffcc00", "Other"="#fdcce5")
thr_pref <- "threshold_quantile"
path_thresholds <- c(glue::glue("{thr_pref}_0.25"), glue::glue("{thr_pref}_0.5"), glue::glue("{thr_pref}_0.75"), glue::glue("{thr_pref}_1"), glue::glue("{thr_pref}_2"),
																					glue::glue("{thr_pref}_3"), glue::glue("{thr_pref}_4"), glue::glue("{thr_pref}_5"), glue::glue("{thr_pref}_10"), glue::glue("{thr_pref}_25"))

# Load metadata and tree
# sc2_md_curated <- readRDS("rds/sc2_md_curated_WITH_Xs_Ns.rds")
# sc2_tre_curated <- readRDS("rds/sc2_tre_curated.rds")
sc2_md_curated <- readRDS(url("https://osf.io/f953r/download", "rb"))
sc2_tre_curated <- readRDS(url("https://osf.io/24r3e/download", "rb"))

#options(max.print=999999)

## 1. Run mlsclust ##
# Period 1: testing (without artifact removal: ~10 minutes: 6 for mlsclust and 4 for run_diff_thresholds)
start <- Sys.time()
res_p1 <- mlsclust(sc2_tre_curated, sc2_md_curated, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, 
																			min_date=as.Date("2019-12-30"), max_date=as.Date("2020-12-31"), branch_length_unit="days", rm_seq_artifacts=TRUE, defining_mut_threshold=0.75, 
																			root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, detailed_output=FALSE, ncpu=NCPU) #detailed_output=TRUE
end <- Sys.time()
total_time_p1 <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p1,"mins"))
if (!dir.exists("rds/")) {
	dir.create("rds/")
}
saveRDS(res_p1, "rds/res_p1.rds")

#res_p1 <- readRDS("rds/res_p1.rds")

start <- Sys.time()
for(i in 1:length(thr)) {
	run_diff_thresholds(stats_df_unfilt=res_p1[[1]], tgt_nodes=res_p1[[2]], homoplasies=res_p1[[3]], 
																					output_dir=glue("results/01_sc2_root_to_dec2020/threshold_quantile_{thr[[i]]}/"), 
																					quantile_choice=quantl[i], quantile_threshold_ratio_sizes=thr[i], quantile_threshold_ratio_persist_time=thr[i], 
																					quantile_threshold_logit_growth=thr[i], threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, 
																					detailed_output=FALSE) #, desc_sisters=list(res_p1[[4]], res_p1[[5]]), amd=res_p1[[6]]
}
end <- Sys.time()
total_time_p1_t <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p1_t,"mins"))

# IMPORTANT: this part is commented because it is time-consuming (see README.md and run in a server if you want to get `res_p2` by yourself)
# Period 2: june 2020 to mid-November 2021
# start <- Sys.time()
# res_p2 <- mlsclust(sc2_tre_curated, sc2_md_curated, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, 
# 																			min_date=as.Date("2019-12-30"), max_date=as.Date("2021-11-15"), branch_length_unit="days", rm_seq_artifacts=TRUE, defining_mut_threshold=0.75, 
# 																			root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, detailed_output=FALSE, ncpu=15)
# end <- Sys.time()
# total_time_p2 <- as.numeric (end - start, units = "mins")
# message(paste("Total time elapsed:",total_time_p2,"mins")) # ~15 hours
#saveRDS(res_p2, "rds/res_p2.rds")
# res_p2 <- readRDS("rds/res_p2.rds") #readRDS(system.file("extdata", "res_p2.rds", package="mlscluster"))
res_p2 <- readRDS(url("https://osf.io/hx9d3/download", "rb"))

start <- Sys.time()
for(i in 1:length(thr)) {
	run_diff_thresholds(stats_df_unfilt=res_p2[[1]], tgt_nodes=res_p2[[2]], homoplasies=res_p2[[3]], 
																					output_dir=glue("results/02_sc2_root_to_nov2021/threshold_quantile_{thr[[i]]}/"), 
																					quantile_choice=quantl[i], quantile_threshold_ratio_sizes=thr[i], quantile_threshold_ratio_persist_time=thr[i], 
																					quantile_threshold_logit_growth=thr[i], threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, 
																					detailed_output=FALSE) #desc_sisters=list(res_p2[[4]], res_p2[[5]]), amd=res_p2[[6]]
}
end <- Sys.time()
total_time_p2_t <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p2_t,"mins"))

# IMPORTANT: this part is commented because it is time-consuming (see README.md and run in a server if you want to get `res_p3` by yourself)
# Period 3: june 2020 to april 2022
# start <- Sys.time()
# res_p3 <- mlsclust(sc2_tre_curated, sc2_md_curated, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, 
# 																			min_date=as.Date("2019-12-30"), max_date=as.Date("2022-04-30"), branch_length_unit="days", rm_seq_artifacts=TRUE, defining_mut_threshold=0.75, 
# 																			root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, detailed_output=FALSE, ncpu=15)
# end <- Sys.time()
# total_time_p3 <- as.numeric (end - start, units = "mins")
# message(paste("Total time elapsed:",total_time_p3,"mins")) # ~44 hours
#saveRDS(res_p3, "rds/res_p3.rds")
# res_p3 <- readRDS("rds/res_p3.rds") #readRDS(system.file("extdata", "res_p3.rds", package="mlscluster"))
res_p3 <- readRDS(url("https://osf.io/9tnck/download", "rb"))

start <- Sys.time()
for(i in 1:length(thr)) {
	run_diff_thresholds(stats_df_unfilt=res_p3[[1]], tgt_nodes=res_p3[[2]], homoplasies=res_p3[[3]], 
																					output_dir=glue("results/03_sc2_whole_period/threshold_quantile_{thr[[i]]}/"), 
																					quantile_choice=quantl[i], quantile_threshold_ratio_sizes=thr[i], quantile_threshold_ratio_persist_time=thr[i], 
																					quantile_threshold_logit_growth=thr[i], threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, 
																					detailed_output=FALSE) #desc_sisters=list(res_p3[[4]], res_p3[[5]]), amd=res_p3[[6]]
}
end <- Sys.time()
total_time_p3_t <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p3_t,"mins"))

## 2. Run statistical tests for multiple thresholds ##
print("=========")
print("PERIOD 2 (PRE-OMICRON)")
print("=========")
# FOR PERIOD 2
# Expected paths containing results for each quantile threshold
table_names_periods <- c(2)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
out_folder <- "period2"
PERIOD_INTEREST <- 2
major_lineages <- c("B.1.177.*","B.1.1.7","AY.4.*",  "AY.* (non-AY.4.*)","Other") #"Beta_B.1.351","Gamma_P.1"

# Change major_lineage labels
change_lineage_labels <- function(path_res) {
	for(i in 1:length(path_thresholds)) {
		print(path_thresholds[i])
		d1 <- utils::read.csv(glue::glue("{path_res}/{path_thresholds[i]}/clustered_all_df.csv"), header=T)
		d1$major_lineage[d1$major_lineage == "Alpha_B.1.1.7"] <- "B.1.1.7"
		d1$major_lineage[d1$major_lineage == "Delta_AY.4.*"] <- "AY.4.*"
		d1$major_lineage[d1$major_lineage == "Delta_other"] <- "AY.* (non-AY.4.*)"
		d1$major_lineage[d1$major_lineage == "EU1_B.1.177"] <- "B.1.177.*"
		d1$major_lineage[d1$major_lineage == "Omicron_BA.1.*"] <- "BA.1.*"
		d1$major_lineage[d1$major_lineage == "Omicron_BA.2.*"] <- "BA.2.*"
		utils::write.csv(d1, file=glue::glue("{path_res}/{path_thresholds[i]}/clustered_all_df.csv"), quote = F, row.names = F)
	}
}

change_lineage_labels("results/02_sc2_root_to_nov2021")

system(glue("mkdir -p stat_results/plots_paper/"))
# Override pal_lineages from package with new major_lineage labels
rmult_p2 <- stats_multiple_thresholds("results/02_sc2_root_to_nov2021", pal_lineages, "period2","ED_FileS1.txt")
plot_tbc_stats("results/02_sc2_root_to_nov2021",out_folder)
stacked_nsites_genomic_region_mult_thresholds("stat_results/period2/genomewide_plot_non-syn/")

print("=========")
print("PERIOD 3 (INCLUDING OMICRON)")
print("=========")
# FOR PERIOD 3
table_names_periods <- c(3)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
out_folder <- "period3"
PERIOD_INTEREST <- 3
#major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Omicron_BA.1.*","Omicron_BA.2.*","Other") #"Beta_B.1.351","Gamma_P.1","Omicron_BA.*",
major_lineages <- c("B.1.177.*","B.1.1.7","AY.4.*", "AY.* (non-AY.4.*)","BA.1.*","BA.2.*","Other")

change_lineage_labels("results/03_sc2_whole_period")

rmult_p3 <- stats_multiple_thresholds("results/03_sc2_whole_period", pal_lineages, "period3","ED_FileS2.txt")

# Extended data Fig S2. freq of syn TFPs along genome 
ed_f2 <- ggarrange(rmult_p3[[3]][[2]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), rmult_p3[[3]][[1]], align='h', nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = FALSE)
ggsave("stat_results/plots_paper/ED_FigS2.png", plot=ed_f2, units="px", width=2200, height=2000, dpi=300, bg="white")

s4 <- plot_tbc_stats("results/03_sc2_whole_period",out_folder)
ggsave("stat_results/plots_paper/ED_FigS4.png", plot=s4, units="px", width=1500, height=2000, dpi=300, bg="white")

# Extended Data Figs S3-S8 (top30 regardless of threshold for main major lineages: AY.4, B.1.1.7, AY.*, BA.1, BA.2)
major_lineages_supp_indices <- c(3,2,4,5,6) #c(5,3,2,4,7,6)
j <- 6
for(i in major_lineages_supp_indices) {
	print(i)
	ggsave(glue("stat_results/plots_paper/ED_FigS{j}_{major_lineages[i]}.png"), plot=rmult_p3[[1]][[i]]+theme(legend.key.size = unit(0.5, 'cm')), units="px", width=2000, height=2200, dpi=300, bg="white")
	j <- j+1
}
# rmult_p3[[1]][[1]] #B.1.177.*
# rmult_p3[[1]][[7]] #Other lineages

sngrt_r3 <- stacked_nsites_genomic_region_mult_thresholds("stat_results/period3/genomewide_plot_non-syn/")

## End run for multiple thresholds ##

## 3. Run FDR ##
# Extracting mutations from the entire dataset of >1.2 million sequences from England
# cog_md_muts_p2 <- extract_muts_period("rds/sc2_md_curated_WITH_Xs_Ns.rds", as.Date("2021-11-15")) # this might need >16 GB of RAM
# cog_md_muts_p3 <- extract_muts_period("rds/sc2_md_curated_WITH_Xs_Ns.rds", as.Date("2022-04-30"))
cog_md_muts_p2 <- extract_muts_period(url("https://osf.io/f953r/download", "rb"), as.Date("2021-11-15")) # this might need >16 GB of RAM
cog_md_muts_p3 <- extract_muts_period(url("https://osf.io/f953r/download", "rb"), as.Date("2022-04-30"))

# Computing codon positions of mutations
muts_match_codons_p2 <- compute_muts_match_codons(cog_md_muts_p2[[1]], cog_md_muts_p2[[2]])
muts_match_codons_p3 <- compute_muts_match_codons(cog_md_muts_p3[[1]], cog_md_muts_p3[[2]])

fdr_plots_p2 <- fdr_mlsclust("results/02_sc2_root_to_nov2021", muts_match_codons_p2[[1]], muts_match_codons_p2[[2]], muts_match_codons_p2[[3]], "period2")
fdr_plots_p3 <- fdr_mlsclust("results/03_sc2_whole_period", muts_match_codons_p3[[1]], muts_match_codons_p3[[2]], muts_match_codons_p3[[3]], "period3")

rem_labels <- function(ggplot_obj_list) {
	for(i in 1:length(ggplot_obj_list)) {
		if(i <= 20) {
			ggplot_obj_list[[i]] <- ggplot_obj_list[[i]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
		}
		if(!(i %in% seq(1,25,5))) {
			ggplot_obj_list[[i]] <- ggplot_obj_list[[i]] + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
		}
		ggplot_obj_list[[i]] <- ggplot_obj_list[[i]] + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			labs(color="Error estimation")
	}
	ggplot_obj_list
}

fdr_plots_p2 <- rem_labels(fdr_plots_p2)
ggarrange(plotlist=fdr_plots_p2, ncol=5, nrow=5, common.legend = T)
ggsave(glue("stat_results/period2/fdr_codons/fdr_all_regions_p2.png"), width=10, height=10, dpi=600, bg="white")

# Extended data Fig S3. FDR and epsilon error rates overall and for each genomic region 
fdr_plots_p3 <- rem_labels(fdr_plots_p3)
fdr_plots_p3_all <- ggarrange(plotlist=fdr_plots_p3, ncol=5, nrow=5, common.legend = T)
annotate_figure(fdr_plots_p3_all, left = text_grob("FDR estimates (%)", rot = 90, vjust = 1, size=8), bottom = text_grob("Thresholds", size=8))
ggsave(glue("stat_results/plots_paper/ED_FigS3.png"), units="px", width=2200, height=2500, dpi=300, bg="white")
## End run FDR ##

## 4. Run for the selected threshold (2% since FDR ~10% and epsilon ~2.5%) ##
print("=========")
print("PERIOD 2 (PRE-OMICRON)")
print("=========")
# FOR PERIOD 2
# Expected paths containing results for each quantile threshold
out_folder <- "period2_thr2"
PERIOD_INTEREST <- 2
THRESHOLD_INTEREST <- 2
#major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Other")

r2 <- stats_selected_threshold("results/02_sc2_root_to_nov2021", thr_index=5, pal_lineages, out_folder) # NOTE: thr_index = 5 -> 2%
# Extended data fig. S5: mlscluster, intersection, and hyphy (also stratified by lineage) identified sites for thr=2% and period before Omicron
r2[[4]]
ggsave(glue("stat_results/plots_paper/ED_FigS5.png"), units="px", width=2200, height=2500, dpi=300, bg="white")


print("=========")
print("PERIOD 3 (INCLUDING OMICRON)")
print("=========")
# FOR PERIOD 3
thr_pref <- "threshold_quantile"
out_folder <- "period3_thr2"
PERIOD_INTEREST <- 3
THRESHOLD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Omicron_BA.1.*","Omicron_BA.2.*","Other")

r3 <- stats_selected_threshold("results/03_sc2_whole_period", thr_index=5, pal_lineages, out_folder)
# Table 1: abs_freq/nonsyn_top30_muts.csv; Table S1: abs_freq/nonsyn_top100_muts.csv

# Figure 3A: number of identified sites for each genomic position stacked by threshold
f5_configs <- theme(axis.title.x=element_text(size=8), axis.title.y=element_text(size=8), axis.text.y=element_text(size=6, color="black"), axis.text.x=element_text(size=7, color="black"))
sngrt_r3 + f5_configs
# Figure 3B: stacked bar counts of unique homoplasies total & norm
leg <- theme(legend.key.size = unit(0.4, 'cm'), legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.75, 'cm'), legend.title = element_text(size=10),legend.text = element_text(size=8)) #guides(color = guide_legend(override.aes = list(size = 1)))
r2[[2]][[1]] <- r2[[2]][[1]] + scale_x_continuous(limits=c(0, 110)) + leg
r3[[2]][[1]] <- r3[[2]][[1]] + scale_x_continuous(limits=c(0, 110)) + leg
rm_titles_all <- theme(axis.title.x=element_blank(), axis.title.y=element_blank())
rm_titles_y <- theme(axis.title.y=element_blank())
rm_axis_text_x <- theme(axis.text.x=element_blank())
my_legend <- get_legend(r3[[2]][[1]])
marg <- ggplot2::theme(plot.margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "cm"), axis.text.y = element_text(size=4))
f3 <- ggarrange(r2[[2]][[1]] + rm_titles_all + rm_axis_text_x + marg, r2[[2]][[2]] + rm_titles_all + rm_axis_text_x + marg, r3[[2]][[1]] + rm_titles_y + marg, r3[[2]][[2]] + rm_titles_y + marg, align='h', vjust=0.5, nrow=2, ncol=2, labels=c('B','C','D','E'), legend.grob=my_legend, legend="right", font.label=list(size=10)) #common.legend = T
annotate_figure(f3, left = text_grob("Genomic region", rot = 90, vjust = 1, size=8))
ggarrange(sngrt_r3 + f5_configs + leg, f3, labels=c("A",""), nrow=2, ncol=1, heights = c(0.4,0.6), widths=c(0.8,1), font.label=list(size=10))
ggsave("stat_results/plots_paper/Fig3.pdf", units="px", width=2100, height=1850, dpi=300, bg="white")

# Figure 4: A. comparison against positive selection (HyPhy) and B. identified sites by major lineage
r3[[4]] #+ theme(axis.text=element_text(size=5))
ggsave("stat_results/plots_paper/Fig4.pdf", units="px", width=2100, height=1850, dpi=300, bg="white")

# Figure 5: spike-wide frequency of TFP-homoplasies
f5_configs <- theme(axis.title.x=element_text(size=8), axis.title.y=element_text(size=7), axis.text.y=element_text(size=6, color="black"), axis.text.x=element_text(size=7, color="black"))
ggarrange(r3[[3]][[10]] + f5_configs, r3[[3]][[3]] + f5_configs, r3[[3]][[6]] + f5_configs, nrow=3, ncol=1, labels=c('A', 'B','C'), legend="right", common.legend = T,  font.label=list(size=10))
ggsave("stat_results/plots_paper/Fig5.pdf", units="px", width=2000, height=1750, dpi=300, bg="white")

# Extended Data Fig. S11: ORF7A and ORF 8 frequencies
ggarrange(r3[[3]][[8]] + f5_configs, r3[[3]][[9]] + f5_configs, nrow=2, ncol=1, labels=c('A', 'B'), legend="right", common.legend = T,  font.label=list(size=10))
ggsave("stat_results/plots_paper/ED_FigS11.png", units="px", width=2000, height=1400, dpi=300, bg="white")

# Extended Data Fig. S1: results for spike when performing 99% quantile frequency outlier removal instead of later adopted alignment-aware method (NOTE: not reproducible here because relying on previous version of the package)
# change_lineage_labels("results/99_quant_rm/03_sc2_whole_period")
# r3_99_rm <- stats_selected_threshold("results/99_quant_rm/03_sc2_whole_period", thr_index=5, pal_lineages, "99_quant_rm/period3_thr2")
# s1_ab <- ggarrange(r3_99_rm[[2]][[1]], r3_99_rm[[2]][[2]], align='h', vjust=0.5, nrow=1, ncol=2, labels=c('A','B'), legend="right", font.label=list(size=10), common.legend = T)
# s1_c <- r3_99_rm[[3]][[10]] + ylim(2,10)
# ggarrange(s1_ab, s1_c, nrow=2, ncol=1, labels=c('A','B'), font.label = list(size=10), common.legend = T)
# ggsave("stat_results/plots_paper/ED_FigS1.png", units="px", width=2000, height=1400, dpi=300, bg="white")
