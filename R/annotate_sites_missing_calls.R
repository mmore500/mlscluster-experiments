# Recovered from old commit here: https://github.com/vinibfranc/mlscluster/blob/cdd6c6ce82b284da31d8be8f8ab9498ab775de8b/R/utils.R
# Found non-syn g10 file here: /home/vinibfranc/Imperial_training/mlscluster-experiments/stat_results/03_BKP_reprod_paper/period2/bef_artifact_removal_old_vers/period3
libs_load <- c("data.table","ape","glue","mlscluster", "readr", "dplyr", "pbmcapply", "stringr")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 3

# Extract potential artifacts for alignment-based analysis (Freq >= 2)
extract_pot_artifacts_g2_freq <- function(df, out_folder, mut_type) {
	df_filt <- df[as.numeric(df$Freq_homopl) >= 2,]
	df_filt <- df_filt[order(df_filt$Freq_homopl, decreasing=T),]
	write.csv(df_filt, glue("stat_results/{out_folder}/{mut_type}_freq_g2.csv"), quote=F, row.names=F)
}

# Gather TFPs for all thresholds
thr <- c(0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
table_names_periods <- c(3)
table_combs <- do.call(paste, c(expand.grid(table_names_periods, thr), sep="_"))
change_ids <- function(df) {
	for(i in 1:nrow(df)) {
		if(df[[i,1]] %in% names(table_combs)) {
			idx <- match(df[[i,1]],names(table_combs))
			df[[i,1]] <- table_combs[idx]
		}
	}
	return(df)
}

# Get all is_clustered==1 homoplasies with Freq >= 2 after running mlsclust without artifact removal (seq_rm_artifacts=FALSE)
tfps_all_thr_bef_artifacts <- function(path_stats) {
	
	if(dir.exists(path_stats))
		path <- path_stats
	else
		stop("Directory does not exist!")
	
	#system(glue("mkdir -p stat_results/{out_folder}/"))
	
	joined_mut_sites_clustered_syn <- joined_mut_sites_clustered_non_syn <- list()
	
	for(i in 1:length(thr)) {
		print("====================")
		print(thr[[i]])
		
		# SYN
		joined_mut_sites_clustered_syn[[i]] <- read.csv(glue("{path_stats}/stat_ready_csvs_only_clustered/df_only_clustered_threshold_quantile_{thr[i]}_SYN.csv"), header=T)
		joined_mut_sites_clustered_non_syn[[i]] <- read.csv(glue("{path_stats}/stat_ready_csvs_only_clustered/df_only_clustered_threshold_quantile_{thr[i]}_NON_SYN.csv"), header=T)
		
	}
	
	names(joined_mut_sites_clustered_syn) <- table_combs
	names(joined_mut_sites_clustered_non_syn) <- table_combs
	
	# Inspection plots
	df_syn_clustered_homopl <- rbindlist(joined_mut_sites_clustered_syn, use.names=T, idcol="period_threshold")
	df_non_syn_clustered_homopl <- rbindlist(joined_mut_sites_clustered_non_syn, use.names=T, idcol="period_threshold")
	
	df_syn_clustered_homopl <- change_ids(df_syn_clustered_homopl)
	df_non_syn_clustered_homopl <- change_ids(df_non_syn_clustered_homopl)
	
	# IMPORTANT: ONLY USE WHEN DOING ARTIFACT REMOVAL
	extract_pot_artifacts_g2_freq(df_syn_clustered_homopl, mut_type="syn", "bef_artifact_removal/period3/")
	extract_pot_artifacts_g2_freq(df_non_syn_clustered_homopl, mut_type="non-syn", "bef_artifact_removal/period3/")
	#list(df_syn_clustered_homopl, df_non_syn_clustered_homopl)
}

tfps_all_thr_bef_artifacts("stat_results/bef_artifact_removal/period3/")

# codon2nucleotide_conversion to handle potential artifacts by checking alignment NNNs (based on https://github.com/theosanderson/Codon2Nucleotide/blob/main/src/App.js)
# Generated the following csvs from code above

nonsyn_p3 <- read.csv("stat_results/bef_artifact_removal/period3/non-syn_freq_g2.csv", header=T)
syn_p3 <- read.csv("stat_results/bef_artifact_removal/period3/syn_freq_g2.csv", header=T)

sc2_md_curated <- readRDS("rds/sc2_md_curated.rds")
system("mkdir -p data/insp/")
write.csv(sc2_md_curated$sequence_name, file="data/insp/incl_seqs.txt", quote=F, sep="\n", row.names=F, col.names=F)
# Important to remove header manually from file above
# Subset 2023 aln to include only sequences included in our study
system("seqtk subseq data/insp/cog_2023-03-31_all_alignment.fa data/insp/incl_seqs.txt > data/insp/incl_seqs.fasta")

codon2nucleotide_conversion <- function(df, non_syn=T) {
	
	# load coordinates
	coords_all <- read.csv(system.file("extdata", "codon2nt_conv_all.csv", package="mlscluster"), header=T) #read.csv("config/codon2nt_conv_all.csv", header=T)
	#coords_nsps <- read.csv("config/codon2nt_conv_nsps.csv", header=T)
	
	# Remove duplicates
	df <- df[!duplicated(df$defining_mut),]
	
	# Separate e.g protein=="S", ref_pos_alt=="E484K" and aa_pos==484
	df$protein <- sub("\\:.*", "", df$defining_mut)
	df$ref_pos_alt <- sub('.*:', "", df$defining_mut)
	df$aa_pos <- parse_number(df$ref_pos_alt)
	
	
	if(non_syn) {
		# join with coord df
		df_join <- df %>% inner_join(coords_all, by="protein")
		df_join$start_nt <- (df_join$start + (df_join$aa_pos - 1) * 3) # NOT DOING THE FOLLOWING ANYMORE: subtracting 1 in the end because COG-UK aln is 29902 nt (not 29903), but not sure this will work for every protein
		df_join$start_nt <- ifelse(df_join$protein == "ORF1AB" & (df_join$start_nt > 13468 & df_join$start_nt < 21557), df_join$start_nt-1, df_join$start_nt)
		#View(df_join)
		df_join$end_nt <- df_join$start_nt + 2
		system("mkdir -p data/insp/non-syn/")
		res_ns <- pbmcapply::pbmclapply(1:nrow(df_join), function(n) {
			system(glue("seqkit subseq -r {df_join[n,21]}:{df_join[n,22]} data/insp/incl_seqs.fasta -o data/insp/non-syn/{df_join[n,7]}:{df_join[n,17]}_nt{df_join[n,21]}-{df_join[n,22]}.fasta")) 
		}, mc.cores = NCPU)
	} else { #SYN
		# identify from which region the mutation come from since only SYNSNP shown
		setDT(df); setDT(coords_all)
		df_join <- df[coords_all, on=c("protein==name_cog","aa_pos>=start","aa_pos<=end"), protein := protein]
		# in this case we will extract only the nt site and not the triplet
		system("mkdir -p data/insp/syn/")
		res_s <- pbmcapply::pbmclapply(1:nrow(df_join), function(n) { 
			system(glue("seqkit subseq -r {df_join[n,17]}:{df_join[n,17]} data/insp/incl_seqs.fasta -o data/insp/syn/{df_join[n,7]}:{df_join[n,17]}.fasta"))
		}, mc.cores = NCPU)
	}
	#View(df_join)
}

#codon2nucleotide_conversion(nonsyn_p2, non_syn=T)
codon2nucleotide_conversion(nonsyn_p3, non_syn=T) #363

#codon2nucleotide_conversion(syn_p2, non_syn=F)
codon2nucleotide_conversion(syn_p3, non_syn=F) #329

dnabin2df <- function(x) {
	df <- data.frame(sequence_name=labels(x),
																		#sp=attr(x, "species"),
																		seq=sapply(as.character(x), paste, collapse=""))
	return(df)
}

# WORKS, BUT MEMORY ISSUES
load_aln_fill_missing_sites <- function(path_ns, path_s) {
	if(dir.exists(path_ns) & dir.exists(path_s)) {
		# extract path for all files
		files_ns <- list.files(path_ns, full.names=T, pattern="*.fasta")

		print("LOADING NON-SYN ALNS")
		
		codon_seqs <- list()
		for(n in 1:length(files_ns)) {
			print(n)
			codon_seqs[[n]] <- read.FASTA(glue("{files_ns[n]}"), type="DNA")
			codon_seqs[[n]] <- dnabin2df(codon_seqs[[n]])
			rownames(codon_seqs[[n]]) <- NULL
			codon_seqs[[n]] <- codon_seqs[[n]][ !(substr(codon_seqs[[n]]$seq,1,1) %in% c("a","c","t","g")) | !(substr(codon_seqs[[n]]$seq,2,2) %in% c("a","c","t","g")) | !(substr(codon_seqs[[n]]$seq,3,3) %in% c("a","c","t","g")),] # if any non-ACTG is found in the codon will add X
			#print(nrow(codon_seqs[[n]]))
		}
		
		# codon_seqs <- pbmcapply::pbmclapply(1:length(files_ns), function(n) { read.FASTA(glue("{files_ns[n]}"), type="DNA") }, mc.cores = NCPU)
		# codon_seqs <- pbmcapply::pbmclapply(1:length(files_ns), function(n) { dnabin2df(codon_seqs[[n]]) }, mc.cores = NCPU)
		# for(n in 1:length(files_ns)) { 
		# 	rownames(codon_seqs[[n]]) <- NULL 
		# 	codon_seqs[[n]] <- codon_seqs[[n]][ !(substr(codon_seqs[[n]]$seq,1,1) %in% c("a","c","t","g")) | !(substr(codon_seqs[[n]]$seq,2,2) %in% c("a","c","t","g")) | !(substr(codon_seqs[[n]]$seq,3,3) %in% c("a","c","t","g")),] # if any non-ACTG is found in the codon will add X
		# }
		
		gc()
		
		# NON-SYN
		print("STARTING NON-SYN MATCHING")
		
		files_basen_df <- dfs_join <- intersect_seqnames <- adj_df1 <- list()
		files_basen <- c()
		for(j in 1:length(files_ns)) {
			print(j)
			# codon_seqs[[j]] <- codon_seqs[[j]][codon_seqs[[j]]$seq == "nnn",]
			#codon_seqs[[j]] <- codon_seqs[[j]][ !(substr(codon_seqs[[j]]$seq,1,1) %in% c("a","c","t","g")) | !(substr(codon_seqs[[j]]$seq,2,2) %in% c("a","c","t","g")) | !(substr(codon_seqs[[j]]$seq,3,3) %in% c("a","c","t","g")),] # if any non-ACTG is found in the codon will add X
			files_basen[j] <- basename(files_ns[j])
			files_basen[j] <- gsub("_.*","",files_basen[j])
			files_basen_df[[j]] <- data.frame(prot_site=files_basen[[j]])
			dfs_join[[j]] <- nonsyn_p3 %>% inner_join(files_basen_df[[j]], by="prot_site")
			dfs_join[[j]]$anc <- str_sub(dfs_join[[j]]$defining_mut, start=1, end=-2)
			intersect_seqnames[[j]] <- sc2_md_curated %>% inner_join(codon_seqs[[j]], by="sequence_name")
			adj_df1[[j]] <- sc2_md_curated %>% inner_join(intersect_seqnames[[j]], by="sequence_name") #left_join
			# !(substr(seq,1,1) %in% c("a","c","t","g")) | !(substr(seq,2,2) %in% c("a","c","t","g")) | !(substr(seq,3,3) %in% c("a","c","t","g")), 
			adj_df1[[j]] <- adj_df1[[j]] %>% mutate(nnn_mutations = if_else(!(substr(seq,1,1) %in% c("a","c","t","g")) | !(substr(seq,2,2) %in% c("a","c","t","g")) | !(substr(seq,3,3) %in% c("a","c","t","g")),
																																																																			paste0("|",unique(dfs_join[[j]]$anc),"X"), mutations.x))
			adj_df1[[j]] <- adj_df1[[j]] %>% select(sequence_name, sample_date.x, lineage.x, major_lineage.x, mutations.x, seq, nnn_mutations) #region.x,
			colnames(adj_df1[[j]]) <- c("sequence_name","sample_date","lineage","major_lineage","mutations","seq","n_mutations") #"region"
		}
		
		rm(codon_seqs, files_basen, intersect_seqnames)
		
		gc()
		
		print("LOADING SYN ALNS")
		
		files_s <- list.files(path_s, full.names=T, pattern="*.fasta")
		
		# nt_seqs <- pbmcapply::pbmclapply(1:length(files_s), function(n) { read.FASTA(glue("{files_s[n]}"), type="DNA") }, mc.cores = NCPU)
		# nt_seqs <- pbmcapply::pbmclapply(1:length(files_s), function(n) { dnabin2df(nt_seqs[[n]]) }, mc.cores = NCPU)
		# for(n in 1:length(files_s)) { 
		# 	rownames(nt_seqs[[n]]) <- NULL 
		# 	nt_seqs[[n]] <- nt_seqs[[n]][!nt_seqs[[n]]$seq %in% c("a","c","t","g"),] #if any non-ACTG is found in base will add N
		# }
		
		nt_seqs <- list()
		for(n in 1:length(files_s)) {
			print(n)
			nt_seqs[[n]] <- read.FASTA(glue("{files_s[n]}"), type="DNA")
			nt_seqs[[n]] <- dnabin2df(nt_seqs[[n]])
			rownames(nt_seqs[[n]]) <- NULL
			nt_seqs[[n]] <- nt_seqs[[n]][!nt_seqs[[n]]$seq %in% c("a","c","t","g"),] #if any non-ACTG is found in base will add N
		}

		gc()
		
		# SYN
		coords_all <- coords_all <- read.csv(system.file("extdata", "codon2nt_conv_all.csv", package="mlscluster"), header = T)
		setDT(syn_p3); setDT(coords_all)
		syn_p3 <- syn_p3[coords_all, on=c("protein==name_cog","site>=start","site<=end"), protein := protein]
		syn_p3$prot_site <- paste0(syn_p3$protein,":",syn_p3$site)
		#syn_p3$defining_mut2 <- gsub("SYNSNP:","synSNP:",syn_p3$defining_mut)
		syn_p3$defining_mut2 <- gsub("SYNSNP:","",syn_p3$defining_mut)
		syn_p3$defining_mut3 <- paste0(syn_p3$protein,":",syn_p3$defining_mut2)
		
		print("STARTING SYN MATCHING")
		
		files_basen2_df <- dfs_join2 <- intersect_seqnames2 <- adj_df1s <- list()
		files_basen2 <- c()
		for(k in 1:length(files_s)) {
			print(k)
			#nt_seqs[[k]] <- nt_seqs[[k]][nt_seqs[[k]]$seq == "n",]
			#nt_seqs[[k]] <- nt_seqs[[k]][!nt_seqs[[k]]$seq %in% c("a","c","t","g"),] #if any non-ACTG is found in base will add N
			files_basen2[k] <- basename(files_s[k])
			files_basen2[k] <- gsub("\\.fasta","",files_basen2[k])
			files_basen2_df[[k]] <- data.frame(prot_site=files_basen2[[k]])
			
			dfs_join2[[k]] <- syn_p3 %>% inner_join(files_basen2_df[[k]], by="prot_site")
			dfs_join2[[k]]$anc <- str_sub(dfs_join2[[k]]$defining_mut3, start=1, end=-2) #defining_mut3
			intersect_seqnames2[[k]] <- sc2_md_curated %>% inner_join(nt_seqs[[k]], by="sequence_name")
			adj_df1s[[k]] <- sc2_md_curated %>% inner_join(intersect_seqnames2[[k]], by="sequence_name") #left_join
			dfs_join2[[k]]$anc <- gsub(".*:","",dfs_join2[[k]]$anc)
			dfs_join2[[k]]$anc <- paste0("synSNP:",dfs_join2[[k]]$anc)
			#seq == 'n', 
			adj_df1s[[k]] <- adj_df1s[[k]] %>% mutate(n_mutations = if_else(!seq %in% c("a","c","t","g"), 
																																																																			paste0("|",unique(dfs_join2[[k]]$anc),"N"), mutations.x))
			adj_df1s[[k]] <- adj_df1s[[k]] %>% select(sequence_name, sample_date.x, lineage.x, major_lineage.x, mutations.x, seq, n_mutations) #region.x,
			colnames(adj_df1s[[k]]) <- c("sequence_name","sample_date","lineage","major_lineage","mutations","seq","n_mutations") #"region"
		}
		
		rm(nt_seqs, files_basen2, dfs_join2, intersect_seqnames2)
		gc()
		
		adj_df2ns <- rbindlist(adj_df1)
		rm(adj_df1); gc()
		adj_df2s <- rbindlist(adj_df1s)
		rm(adj_df1s); gc()
		adj_df2 <- rbind(adj_df2ns, adj_df2s)
		rm(adj_df2ns, adj_df2s); gc()
		adj_df2 <- adj_df2 %>% group_by(sequence_name) %>% summarize(miss_mutations = paste0(na.omit(n_mutations), collapse="")) %>% ungroup()
		adj_df3 <- sc2_md_curated %>% left_join(adj_df2, by="sequence_name")
		rm(adj_df2); gc()
		adj_df3$muts_joined <- paste0(adj_df3$mutations, adj_df3$miss_mutations) #na.omit
		adj_df3 <- adj_df3 %>% select(sequence_name, sample_date, lineage, major_lineage, muts_joined) #region
		colnames(adj_df3) <- c("sequence_name", "sample_date", "lineage", "major_lineage", "mutations") #"region",
		adj_df3$mutations <- gsub("NA","", as.character(adj_df3$mutations))
	}else {
		print("One or both the specified paths do not exist!")
	}
	return(adj_df3)
}

incl_missing_muts <- load_aln_fill_missing_sites("data/insp/non-syn", "data/insp/syn")
incl_missing_muts$mutations <- toupper(incl_missing_muts$mutations)
#sc2_md_curated <- incl_missing_muts
saveRDS(incl_missing_muts, file="rds/sc2_md_curated_WITH_Xs_Ns.rds")

# incl_missing_muts[incl_missing_muts$sequence_name == "England/ALDP-13A2064/2021",]$mutations #both syn and non-syn
# incl_missing_muts[incl_missing_muts$sequence_name == "England/ALDP-14A91A5/2021",]$mutations #only syn
# incl_missing_muts[incl_missing_muts$sequence_name == "England/ALDP-14ADC37/2021",]$mutations #only non-syn