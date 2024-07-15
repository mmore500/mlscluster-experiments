# Move DMS files from https://github.com/jbloomlab/SARS2-mut-fitness/tree/main/results/dms to data/dms
# Move FITNESS files from https://github.com/jbloomlab/SARS2-mut-fitness/tree/main/results_public_2022-01-31 to data/fitness

libs_load <- c("glue","dplyr", "ggplot2","ggforce", "ggrepel", "ggpubr")
invisible( lapply(libs_load, library, character.only=TRUE) )

### Compare against deep-mutational scanning (DMS) experimental results 

dms_files <- list.files(path="data/dms", full.names = T)

# Adjust DMS files

dms_csv <- list()
for(i in 1:length(dms_files)) {
	dms_csv[[i]] <- read.csv(dms_files[i], header=T)
	dms_csv[[i]]$site <- as.integer(dms_csv[[i]]$site)
	#print(head(dms_csv))
	if (grepl("spike", dms_files[i]) | grepl("rbd", dms_files[i])) {
		dms_csv[[i]]$protein <- "S"
		print("1")
		print(dms_files[i])
	} else if(grepl("mpro", dms_files[i])) {
		dms_csv[[i]]$protein <- "NSP5"
		print("2")
		print(dms_files[i])
	}
	dms_csv[[i]]$defining_mut <- glue("{dms_csv[[i]]$protein}:{dms_csv[[i]]$wildtype}{dms_csv[[i]]$site}{dms_csv[[i]]$mutant}")
	print(range(dms_csv[[i]]$site, na.rm=T))
	#View(dms_csv)
}
#dms_csv[[1]]$site 214a, b and c???

dms_csv[[5]]$s_mut_region_interest <- "RBD"

# Read tfps
folder_all_tfps_t2 <- "results/03_sc2_whole_period/threshold_quantile_2"
# Load tfps and non-tfps
tfps_t2_all <- read.csv(glue("{folder_all_tfps_t2}/clustered_all_df.csv"), header=T)
tfps_t2 <- tfps_t2_all[tfps_t2_all$is_clustered == 1,]
# non-syn tfps
tfps_t2_ns <- tfps_t2[tfps_t2$protein != "SYNSNP",]
# Adjust NSP5 coordinates
tfps_t2_ns$site <- sub('.*:', "", tfps_t2_ns$defining_mut)
tfps_t2_ns$site <- readr::parse_number(tfps_t2_ns$site, na="X")
tfps_t2_ns$site <- as.integer(tfps_t2_ns$site)
tfps_t2_ns$anc_site_mut <- sub('.*:', "", tfps_t2_ns$defining_mut)
tfps_t2_ns$mutsite <- stringr::str_sub(tfps_t2_ns$defining_mut, -1)
tfps_t2_ns$ancsite <- stringr::str_sub(tfps_t2_ns$anc_site_mut, start=1, end=1)

tfps_t2_ns$defining_mut <- ifelse( (tfps_t2_ns$protein == "ORF1AB") & (tfps_t2_ns$site >= 3264 & tfps_t2_ns$site <= 3570), 
																																			paste0("NSP5:",tfps_t2_ns$ancsite, "", tfps_t2_ns$site - 3263, "", tfps_t2_ns$mutsite), tfps_t2_ns$defining_mut)
tfps_t2_ns$protein <- ifelse( (tfps_t2_ns$protein == "ORF1AB") & (tfps_t2_ns$site >= 3264 & tfps_t2_ns$site <= 3570), "NSP5", tfps_t2_ns$protein)

# non-syn non-tfps
non_tfps_t2 <- tfps_t2_all[tfps_t2_all$is_clustered == 0,]
non_tfps_t2_ns <- non_tfps_t2[non_tfps_t2$protein != "SYNSNP",]
non_tfps_t2_ns$site <- sub('.*:', "", non_tfps_t2_ns$defining_mut)
non_tfps_t2_ns$site <- readr::parse_number(non_tfps_t2_ns$site, na="X")
non_tfps_t2_ns$site <- as.integer(non_tfps_t2_ns$site)
non_tfps_t2_ns$anc_site_mut <- sub('.*:', "", non_tfps_t2_ns$defining_mut)
non_tfps_t2_ns$mutsite <- stringr::str_sub(non_tfps_t2_ns$defining_mut, -1)
non_tfps_t2_ns$ancsite <- stringr::str_sub(non_tfps_t2_ns$anc_site_mut, start=1, end=1)

non_tfps_t2_ns$defining_mut <- ifelse( (non_tfps_t2_ns$protein == "ORF1AB") & (non_tfps_t2_ns$site >= 3264 & non_tfps_t2_ns$site <= 3570), 
																																			paste0("NSP5:",non_tfps_t2_ns$ancsite, "", non_tfps_t2_ns$site - 3263, "", non_tfps_t2_ns$mutsite), non_tfps_t2_ns$defining_mut)
non_tfps_t2_ns$protein <- ifelse( (non_tfps_t2_ns$protein == "ORF1AB") & (non_tfps_t2_ns$site >= 3264 & non_tfps_t2_ns$site <= 3570), "NSP5", non_tfps_t2_ns$protein)


match_dms_tfps <- function(dms_meas, tfps, non_tfps, prot_interest, rbd_only=F) {

	print("Distribution of both positive and negative effects for all DMS data: ")
	hist(dms_meas$effect, breaks=50)	
	
	print("Total sites with DMS data: ")
	print(nrow(dms_meas))
	dms_meas_neg <- dms_meas[dms_meas$effect < 0,]
	print("Sites with DMS negative effect: ")
	print(nrow(dms_meas_neg))
	print("============")
	dms_meas_pos <- dms_meas[dms_meas$effect > 0,]
	print("Sites with DMS positive effect: ")
	print(nrow(dms_meas_pos))
	print("============")
	dms_meas_neg_prot <- dms_meas_neg[dms_meas_neg$protein == prot_interest,]
	dms_meas_pos_prot <- dms_meas_pos[dms_meas_pos$protein == prot_interest,]
	
	if(rbd_only)  {
		dms_meas_neg_prot <- dms_meas_neg_prot[dms_meas_neg_prot$s_mut_region_interest == "RBD",]
		dms_meas_pos_prot <- dms_meas_pos_prot[dms_meas_pos_prot$s_mut_region_interest == "RBD",]
	}
	
	# print(glue("Sites with negative effects on {prot_interest}"))
	# print(nrow(dms_meas_neg_prot))
	tfps_prot <- tfps[tfps$protein == prot_interest,]
	non_tfps_prot <- non_tfps[non_tfps$protein == prot_interest,]
	
	if(rbd_only) {
		tfps_prot <- tfps_prot[tfps_prot$s_mut_region_interest == "RBD" & !is.na(tfps_prot$s_mut_region_interest),]
		non_tfps_prot <- non_tfps_prot[non_tfps_prot$s_mut_region_interest == "RBD" & !is.na(non_tfps_prot$s_mut_region_interest),]
	}
	print(glue("TFPs on {prot_interest}"))
	print(nrow(tfps_prot))
	
	print(glue("Non-TFPs on {prot_interest}"))
	print(nrow(non_tfps_prot))
	
	matched_tfps_pos <- inner_join(dms_meas_pos, tfps_prot, by="defining_mut")
	print("Number of matching POSITIVE dms estimates and TFP sites: ")
	print(nrow(matched_tfps_pos))
	#print(matched_tfps_pos %>% group_by(protein) %>% summarise(n=n()))
	print("Distribution of POSITIVE effects for matched fitness estimates and TFPs: ")
	hist(matched_tfps_pos$effect, breaks=50)
	
	matched_non_tfps_pos <- inner_join(dms_meas_pos, non_tfps_prot, by="defining_mut")
	print("Number of matching POSITIVE fitness estimates and NOT TFP sites: ")
	print(nrow(matched_non_tfps_pos))
	#print(matched_non_tfps_pos %>% group_by(protein) %>% summarise(n=n()))
	
	matched_tfps_neg <- inner_join(dms_meas_neg, tfps_prot, by="defining_mut")
	print("Number of matching NEGATIVE dms estimates and TFP sites: ")
	print(nrow(matched_tfps_neg))
	#print(matched_non_tfps_neg %>% group_by(protein) %>% summarise(n=n()))
	
	matched_non_tfps_neg <- inner_join(dms_meas_neg, non_tfps_prot, by="defining_mut")
	print("Number of matching NEGATIVE fitness estimates and NOT TFP sites: ")
	print(nrow(matched_non_tfps_neg))
	#print(matched_non_tfps_neg %>% group_by(protein) %>% summarise(n=n()))
	
	# crosstab (2x2)
	# 										TFP 	Not TFP
	# DMS+	
	# DMS-
	
	crosstab <- matrix(nrow=2, ncol=2)
	rownames(crosstab) <- c("DMS+", "DMS-")
	colnames(crosstab) <- c("TFP", "Not TFP")
	crosstab[[1,1]] <- nrow(matched_tfps_pos)
	crosstab[[2,1]] <- nrow(matched_tfps_neg)
	crosstab[[1,2]] <- nrow(matched_non_tfps_pos)
	crosstab[[2,2]] <- nrow(matched_non_tfps_neg)
	
	print(crosstab)
	
	list(crosstab=crosstab, matched_tfps_pos=matched_tfps_pos, matched_tfps_neg=matched_tfps_neg, matched_non_tfps_pos=matched_non_tfps_pos, matched_non_tfps_neg=matched_non_tfps_neg)
}

# Compare differences of DMS or Bloom selection coefficients between TFPs and non-TFPs
compare_distr_coeff <- function(prev_df, variable, variable_label) {
	#ggplot(data=prev_df, aes(x=as.factor(is_clustered), y=!!sym(variable), color=as.factor(is_clustered))) + scale_color_discrete(name="Is TFP?") + geom_violin() + geom_sina() + labs(x="Is TFP?", y=variable_label) + theme_minimal()
	prev_df$Status[prev_df$is_clustered == 1] <- "TFP"
	prev_df$Status[prev_df$is_clustered == 0] <- "non-TFP"
	pl <- ggplot(data=prev_df, aes(x=Status, y=!!sym(variable))) + geom_violin(color="#3d5a80") + geom_sina(color="#3d5a80") + labs(x="Mutation status", y=variable_label) + theme_minimal() + theme(legend.position = "none", axis.text.x=element_text(size=6, color="black"), axis.text.y=element_text(size=6,color="black"), axis.title=element_text(size=7)) # color=Status, scale_color_discrete(name="Is TFP?")
	pl
}


extract_directionality_diff <- function(df, variable) {
	median_tfp <- median(df[[variable]][df$is_clustered == 1])
	print(paste0("Median ",variable," for TFP: ",median_tfp))
	median_non_tfp <- median(df[[variable]][df$is_clustered == 0])
	print(paste0("Median ",variable," for non-TFP: ",median_non_tfp))
	
	if (median_tfp > median_non_tfp) {
		directionality <- paste0("TFP has higher median ",variable)
	} else if (median_tfp < median_non_tfp) {
		directionality <- paste0("non-TFP has higher median ",variable)
	} else {
		directionality <- paste0("TFP and non-TFP have the same ",variable)
	}
	
	directionality
}

# Dadonaite spike
dms_tfp_spike1 <- match_dms_tfps(dms_csv[[1]], tfps_t2_ns, non_tfps_t2_ns, "S") # 73 TFPs, 49 with negative effects (-2.5 to 0), 15 with slightly positive (max: 0.3483)
fisher.test(dms_tfp_spike1$crosstab)
fisher.test(matrix( c( 49, 15 , 332, 169), ncol=2)) # inverse
hist(dms_tfp_spike1$matched_tfps_neg$effect); range(dms_tfp_spike1$matched_tfps_neg$effect)
dms_tfp_spike1$matched_tfps_neg$defining_mut
hist(dms_tfp_spike1$matched_tfps_pos$effect); range(dms_tfp_spike1$matched_tfps_pos$effect)
dms_tfp_spike1$matched_tfps_pos$defining_mut

pos_dms_s_tfps <- dms_tfp_spike1$matched_tfps_pos %>% select(defining_mut, effect, is_clustered)
pos_dms_s_non_tfps <- dms_tfp_spike1$matched_non_tfps_pos %>% select(defining_mut, effect, is_clustered)
pos_dms_comb <- rbind(pos_dms_s_tfps, pos_dms_s_non_tfps)
compare_distr_coeff(pos_dms_comb, "effect", "DMS effect")
# test diff (signif. p-value indicates significant difference in the distributions of effect between TFP and non-TFP)
wilcox.test(effect ~ is_clustered, data = pos_dms_comb) #p-value = 0.02164
extract_directionality_diff(pos_dms_comb, "effect")

# Starr RBD
dms_tfp_rbd <- match_dms_tfps(dms_csv[[5]], tfps_t2_ns, non_tfps_t2_ns, "S", rbd_only = T) # 5 TFPs, 4 with negative effects, 1 with positive effect (0.078808)
fisher.test(dms_tfp_rbd$crosstab)
fisher.test(matrix( c( 4, 1 , 59, 16), ncol=2)) # inverse
hist(dms_tfp_rbd$matched_tfps_neg$effect); range(dms_tfp_rbd$matched_tfps_neg$effect)
hist(dms_tfp_rbd$matched_tfps_pos$effect); range(dms_tfp_rbd$matched_tfps_pos$effect)

plot_dms_tfp_matches_over_gene <- function(df_overlaps, label, break_intervals_x) {
	# join all to plot
	tfp_dms_pos <- df_overlaps$matched_tfps_pos; tfp_dms_pos$Overlap <- "DMS+ / TFP"
	tfp_dms_neg <- df_overlaps$matched_tfps_neg; tfp_dms_neg$Overlap <- "DMS- / TFP"
	non_tfp_dms_pos <- df_overlaps$matched_non_tfps_pos; non_tfp_dms_pos$Overlap <- "DMS+ / non-TFP"
	non_tfp_dms_neg <- df_overlaps$matched_non_tfps_neg; non_tfp_dms_neg$Overlap <- "DMS- / non-TFP"
	crosstab_df <- rbind(tfp_dms_pos, tfp_dms_neg, non_tfp_dms_pos, non_tfp_dms_neg)
	
	ggplot(crosstab_df, aes(y=effect, x=site.x, color=Overlap)) + #color=drugs
		geom_point(alpha=0.8) + labs(y="DMS effect",x=glue("{label} site")) +
		geom_text_repel(aes(label = ifelse(Overlap %in% c("DMS+ / TFP", "DMS- / TFP"),as.character(anc_site_mut),'')), size=2, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=30) + 
		scale_color_brewer(palette = "Spectral") + scale_x_continuous(breaks = seq(0, max(crosstab_df$site.x), by = break_intervals_x)) +
		theme_minimal() + theme(legend.position="top", axis.text.x=element_text(size=6, color="black"), axis.text.y=element_text(size=6,color="black"), axis.title=element_text(size=7)) #angle=90, vjust = 0.5, hjust=1,
}

pl_s <- plot_dms_tfp_matches_over_gene(dms_tfp_spike1, "Spike", 200) 
pl_rbd <- plot_dms_tfp_matches_over_gene(dms_tfp_rbd, "Spike (RBD)", 25)

ggarrange(pl_s, pl_rbd, nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = TRUE)
# Fig. S12
ggsave(glue("stat_results/plots_paper/ED_FigS12.png"), units="px", width=2200, height=2500, dpi=300, bg="white")

# Flynn mpro (nsp5) 2022
dms_tfp_nsp5_1 <- match_dms_tfps(dms_csv[[2]], tfps_t2_ns, non_tfps_t2_ns, "NSP5") # 12 TFPs, 0 with negative effects: 0.72352 to 1.02030 effect (only 748 out of 6048 sites with negative effects???)
dms_tfp_nsp5_1$matched_tfps_pos$defining_mut
hist(dms_tfp_nsp5_1$matched_tfps_pos$effect); range(dms_tfp_nsp5_1$matched_tfps_pos$effect)
# Flynn mpro (nsp5) 2023
dms_tfp_nsp5_2 <- match_dms_tfps(dms_csv[[3]], tfps_t2_ns, non_tfps_t2_ns, "NSP5") # 12 TFPs, 0 with negative effects (all unmatched): 0.4667 to 0.7350 effect (only 565 out of 6025 sites with negative effects???)
hist(dms_tfp_nsp5_2$matched_tfps_pos$effect); range(dms_tfp_nsp5_2$matched_tfps_pos$effect)
# Iketani mpro (nsp5)
dms_tfp_nsp5_3 <- match_dms_tfps(dms_csv[[4]], tfps_t2_ns, non_tfps_t2_ns, "NSP5") # 12 TFPs, 2 with negative effects, 8 with positive effects: 0.025325 to 0.672010 (3997 out of 6060 sites with negative effects)
hist(dms_tfp_nsp5_3$matched_tfps_pos$effect); range(dms_tfp_nsp5_3$matched_tfps_pos$effect)
hist(dms_tfp_nsp5_3$matched_tfps_neg$effect); range(dms_tfp_nsp5_3$matched_tfps_neg$effect)

### Compare against Bloom & Neher paper

# only England
aa_england_fit_jan2022 <- read.csv("data/fitness/aamut_fitness_by_subset_2022-01-31.txt")
aa_england_fit_jan2022 <- aa_england_fit_jan2022[aa_england_fit_jan2022$subset == "England",]
aa_england_fit_jan2022$defining_mut <- paste0(aa_england_fit_jan2022$gene, ":", aa_england_fit_jan2022$aa_mutation)
aa_england_fit_jan2022$defining_mut <-	toupper(aa_england_fit_jan2022$defining_mut)
aa_england_fit_jan2022 <- aa_england_fit_jan2022[aa_england_fit_jan2022$clade_founder_aa != aa_england_fit_jan2022$mutant_aa,] #121419 to 105117

# Load tfps (same as before but with more annotations)
folder2_all_tfps_t2 <- "stat_results/period3_thr2/stat_ready_csvs_only_clustered/"
tfps_t2_ns_comp <- read.csv(glue("{folder2_all_tfps_t2}/df_only_clustered_threshold_quantile_2_NON_SYN.csv"))

match_bloomNeher_tfps <- function(fit_bloomNeher, tfps, non_tfps, non_tfp_comparison="mlscluster_non_tfps") {
	
	#print("Distribution of effects for both + and - fitness estimates: ")
	#hist(fit_bloomNeher$delta_fitness, breaks=100)
	
	print("Total sites with fitness estimates: ")
	print(nrow(fit_bloomNeher))
	fit_bloomNeher_neg <- fit_bloomNeher[fit_bloomNeher$delta_fitness < 0,]
	print("Sites with negative fitness effect: ")
	print(nrow(fit_bloomNeher_neg))
	fit_bloomNeher_pos <- fit_bloomNeher[fit_bloomNeher$delta_fitness > 0,]
	print("Sites with positive fitness effect: ")
	print(nrow(fit_bloomNeher_pos))
	print("============")
	
	#print("Distribution of effects for all negative fitness estimates: ")
	#hist(fit_bloomNeher_neg$delta_fitness, breaks=50)
	
	print(glue("TFPs: "))
	print(nrow(tfps))
	print(tfps %>% group_by(protein) %>% summarise(n=n()))
	
	matched_tfps_neg <- inner_join(fit_bloomNeher_neg, tfps, by="defining_mut")
	print("Number of matching NEGATIVE fitness estimates and TFP sites: ")
	print(nrow(matched_tfps_neg))
	print(matched_tfps_neg %>% group_by(protein) %>% summarise(n=n()))
	
	#print("Distribution of NEGATIVE effects for matched fitness estimates and TFPs: ")
	#hist(matched_tfps_neg$delta_fitness, breaks=50)
	
	matched_tfps_pos <- inner_join(fit_bloomNeher_pos, tfps, by="defining_mut")
	print("Number of matching POSITIVE fitness estimates and TFP sites: ")
	print(nrow(matched_tfps_pos))
	print(matched_tfps_pos %>% group_by(protein) %>% summarise(n=n()))
	print("Distribution of POSITIVE effects for matched fitness estimates and TFPs: ")
	hist(matched_tfps_pos$delta_fitness, breaks=50)
	
	if(non_tfp_comparison == "mlscluster_non_tfps") {
		matched_non_tfps_neg <- inner_join(fit_bloomNeher_neg, non_tfps, by="defining_mut")
		print("Number of matching NEGATIVE fitness estimates and NOT TFP sites: ")
		print(nrow(matched_non_tfps_neg))
		print(matched_non_tfps_neg %>% group_by(protein) %>% summarise(n=n()))
		
		matched_non_tfps_pos <- inner_join(fit_bloomNeher_pos, non_tfps, by="defining_mut")
		print("Number of matching POSITIVE fitness estimates and NOT TFP sites: ")
		print(nrow(matched_non_tfps_pos))
		print(matched_non_tfps_pos %>% group_by(protein) %>% summarise(n=n()))
	} else {
		matched_non_tfps_neg <- fit_bloomNeher_neg[!(fit_bloomNeher_neg$defining_mut %in% matched_tfps_neg$defining_mut),]
		print("Number of matching NEGATIVE fitness estimates and NOT TFP sites: ")
		print(nrow(matched_non_tfps_neg))
		#print(matched_non_tfps_neg %>% group_by(protein) %>% summarise(n=n()))
		
		matched_non_tfps_pos <- fit_bloomNeher_pos[!(fit_bloomNeher_pos$defining_mut %in% matched_tfps_pos$defining_mut),]
		print("Number of matching POSITIVE fitness estimates and NOT TFP sites: ")
		print(nrow(matched_non_tfps_pos))
		#print(matched_non_tfps_pos %>% group_by(protein) %>% summarise(n=n()))
	}
	
	# crosstab (2x2)
	# 										TFP 	Not TFP
	# Bloom+	
	# Bloom-
	
	crosstab <- matrix(nrow=2, ncol=2)
	rownames(crosstab) <- c("Bloom+", "Bloom-")
	colnames(crosstab) <- c("TFP", "Not TFP")
	crosstab[[1,1]] <- nrow(matched_tfps_pos)
	crosstab[[2,1]] <- nrow(matched_tfps_neg)
	crosstab[[1,2]] <- nrow(matched_non_tfps_pos)
	crosstab[[2,2]] <- nrow(matched_non_tfps_neg)
	
	#return(list(matched_tfps_neg=matched_tfps_neg ))
	print(crosstab)
	
	list(crosstab=crosstab, matched_tfps_pos=matched_tfps_pos, matched_tfps_neg=matched_tfps_neg, matched_non_tfps_pos=matched_non_tfps_pos, matched_non_tfps_neg=matched_non_tfps_neg)
}

# Make sure there are no overlaps between tfp and non-tfp
occur <- table(factor(non_tfps_t2_ns$defining_mut, levels = tfps_t2_ns_comp$defining_mut))
sum(non_tfps_t2_ns$defining_mut %in% tfps_t2_ns_comp$defining_mut)
sum(tfps_t2_ns_comp$defining_mut %in% non_tfps_t2_ns$defining_mut)

matches_england <- match_bloomNeher_tfps(aa_england_fit_jan2022, tfps_t2_ns_comp, non_tfps_t2_ns, non_tfp_comparison = "mlscluster_non_tfps") # of the 547 TFPs, 459 (83.9%) also detected as deleterious by Bloom & Neher for England subset, 80 TFPs have fitness usually between 0 and 1
View(matches_england$matched_tfps_pos)
nrow(matches_england$matched_tfps_pos) #80
nrow(matches_england$matched_tfps_pos[matches_england$matched_tfps_pos$delta_fitness > 2,]) #>1: 10, >2: 4
hist(matches_england$matched_tfps_pos$delta_fitness, breaks=50)
View(matches_england$matched_tfps_neg)
hist(matches_england$matched_tfps_neg$delta_fitness, breaks=50)

fisher.test(matches_england$crosstab) ##about 1.8 higher odds that a TFP will be s>0
fisher.test(matrix(c(459,80,3633,351),ncol=2)) # inverse

chisq_england <- chisq.test(matches_england$crosstab)
print(chisq_england)
print(chisq_england$observed)
print(chisq_england$expected)

# pos
pos_tfps <- matches_england$matched_tfps_pos %>% select(defining_mut, delta_fitness, is_clustered, gene)
pos_non_tfps <- matches_england$matched_non_tfps_pos %>% select(defining_mut, delta_fitness, is_clustered, gene)
matches_england_pos <- rbind(pos_tfps, pos_non_tfps)

pl_fit_pos <- compare_distr_coeff(matches_england_pos, "delta_fitness", "Delta fitness")
hist(matches_england_pos$delta_fitness[matches_england_pos$is_clustered == 1], breaks=100)
hist(matches_england_pos$delta_fitness[matches_england_pos$is_clustered == 0], breaks=100)
# test diff (signif. p-value indicates significant difference in the distributions of selection [delta fitness] between TFP and non-TFP)
wilcox.test(delta_fitness ~ is_clustered, data = matches_england_pos) #p-value = 0.008036

# neg
neg_tfps <- matches_england$matched_tfps_neg %>% select(defining_mut, delta_fitness, is_clustered, gene)
neg_non_tfps <- matches_england$matched_non_tfps_neg %>% select(defining_mut, delta_fitness, is_clustered, gene)
matches_england_neg <- rbind(neg_tfps, neg_non_tfps)

pl_fit_neg <- compare_distr_coeff(matches_england_neg, "delta_fitness", "Delta fitness")
hist(matches_england_neg$delta_fitness[matches_england_neg$is_clustered == 1], breaks=100)
hist(matches_england_neg$delta_fitness[matches_england_neg$is_clustered == 0], breaks=100)
wilcox.test(delta_fitness ~ is_clustered, data = matches_england_neg) #p-value = 0.000000000005307

extract_directionality_diff(matches_england_pos, "delta_fitness")
extract_directionality_diff(matches_england_neg, "delta_fitness")

plot_fitness_bloom_tfp_matches_over_gene <- function(df_overlaps, selected_prot, label, break_intervals_x) {
	# join all to plot
	tfp_fit_pos <- df_overlaps$matched_tfps_pos; tfp_fit_pos$Overlap <- "Fitness+ / TFP"; tfp_fit_pos <- tfp_fit_pos %>% select(defining_mut, site, protein, delta_fitness, Overlap) #anc_site_mut,spec_regions
	tfp_fit_neg <- df_overlaps$matched_tfps_neg; tfp_fit_neg$Overlap <- "Fitness- / TFP"; tfp_fit_neg <- tfp_fit_neg %>% select(defining_mut, site, protein, delta_fitness, Overlap)
	non_tfp_fit_pos <- df_overlaps$matched_non_tfps_pos; non_tfp_fit_pos$Overlap <- "Fitness+ / non-TFP"; non_tfp_fit_pos <- non_tfp_fit_pos %>% select(defining_mut, site, protein, delta_fitness, Overlap)
	non_tfp_fit_neg <- df_overlaps$matched_non_tfps_neg; non_tfp_fit_neg$Overlap <- "Fitness- / non-TFP"; non_tfp_fit_neg <- non_tfp_fit_neg %>% select(defining_mut, site, protein, delta_fitness, Overlap)
	crosstab_df <- rbind(tfp_fit_pos, tfp_fit_neg, non_tfp_fit_pos, non_tfp_fit_neg)
	crosstab_df <- crosstab_df[crosstab_df$protein == selected_prot,]
	crosstab_df$anc_site_mut <- sub('.*:', "", crosstab_df$defining_mut)
	
	pl <- ggplot(crosstab_df, aes(y=delta_fitness, x=site, color=Overlap)) + #color=drugs
		geom_point(alpha=0.8) + labs(y="Delta fitness",x=glue("{label} site")) +
		geom_text_repel(aes(label = ifelse(Overlap %in% c("Fitness+ / TFP", "Fitness- / TFP"),as.character(anc_site_mut),'')), size=1.5, box.padding = 0.2, point.padding = 0.2, segment.color = 'grey50', max.overlaps=30) + #2,0.35,0.5
		scale_color_brewer(palette = "Spectral") + scale_x_continuous(breaks = seq(0, max(crosstab_df$site), by = break_intervals_x)) +
		theme_minimal() + theme(legend.position="top", axis.text.x=element_text(size=6, color="black"), axis.text.y=element_text(size=6,color="black"), axis.title=element_text(size=7)) #angle=90, vjust = 0.5, hjust=1,
	pl
}

pl_fit_s <- plot_fitness_bloom_tfp_matches_over_gene(matches_england, "S", "Spike", break_intervals_x = 200)
pl_fit_n <- plot_fitness_bloom_tfp_matches_over_gene(matches_england, "N", "Nucleocapsid", break_intervals_x = 50)
pl_fit_orf3a <- plot_fitness_bloom_tfp_matches_over_gene(matches_england, "ORF3A", "ORF3a", break_intervals_x = 50)

s13_1 <- ggarrange(pl_fit_pos, pl_fit_neg, ncol=2, nrow=1, labels=c("A","B"))
s13_2 <- ggarrange(pl_fit_s, pl_fit_n, pl_fit_orf3a, ncol=1, nrow=3, labels=c("C","D","E"), common.legend = T)
ggarrange(s13_1, s13_2, heights = c(2,8), ncol=1, nrow=2) #ncol=1, nrow=2

# Fig. S13
#ggsave(glue("stat_results/plots_paper/ED_FigS13.png"), units="px", width=2500, height=1500, dpi=300, bg="white")
ggsave(glue("stat_results/plots_paper/ED_FigS13.png"), units="px", width=2200, height=2700, dpi=300, bg="white")
