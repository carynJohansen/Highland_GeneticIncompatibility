# Sleuth analysis for kallisto output

# basic code found at: https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

setwd("~/Box Sync/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest")

library(sleuth)

sample_id <- dir(file.path("data/fullTranscript/"))
sample_id

#get the file path for the results
kal_results <- file.path("data/fullTranscript", sample_id)
kal_results

res <- as.data.frame(cbind(sample_id, kal_results))

#generate the data structure to describe the samples
meta <- read.table("simple_meta.txt", header=T, sep = " ", stringsAsFactors = T)
meta_leaf <- dplyr::filter(meta, tissue_s == "leaf")
meta_leaf <- dplyr::select(meta_leaf, sample=Run_s, cultivar=cultivar_s)
str(meta_leaf)

sample_leaf <- as.character(meta_leaf$sample)
res_leaf <- res[as.character(res$sample_id) %in% sample_leaf,]
meta_leaf <- dplyr::mutate(meta_leaf, path=as.character(res_leaf$kal_results))

# Generate the sleuth object:

leaf_matrix <- model.matrix(~cultivar, meta_leaf)

leaf_so <- sleuth_prep(meta_leaf, ~cultivar, extra_bootstrap_summary = TRUE)

# full model fit
# to smooth the raw kallisto abundance estimates, using a linear model
# with parameter that represents the experimental conditions
leaf_so <- sleuth_fit(leaf_so, leaf_matrix, 'full')

# this second fit to a reduced model presumes the abundances are equal
# between the two conditions
leaf_so <- sleuth_fit(leaf_so, ~1, 'reduced')

#test
leaf_so <- sleuth_lrt(leaf_so, 'reduced', 'full')

models(leaf_so)

#saveRDS(leaf_so, file="leaf_sleuth_results.rds")

#results of test
sleuth_table <- sleuth_results(leaf_so, 'reduced:full', 'lrt', show_all=FALSE)

#saveRDS(sleuth_table, file="leaf_table.rds")

sleuth_sig <- dplyr::filter(sleuth_table, qval<=0.05)
dim(sleuth_sig)
head(sleuth_sig, 20)


plot_bootstrap(leaf_so, "Zm00001d005814_T001", units="est_counts", color_by="cultivar")

plot_bootstrap(leaf_so, sleuth_sig[2,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[3,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[4,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[5,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[6,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[7,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[8,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[9,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[10,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[11,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[12,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[2000,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[3000,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[5000,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[4900,]$target_id, units="est_counts", color_by="cultivar")
plot_bootstrap(leaf_so, sleuth_sig[2500,]$target_id, units="est_counts", color_by="cultivar")

plot_pca(leaf_so, color_by = 'cultivar')

#Aggregate p-values for gene-based analysis

# as done in 'Gene Level Differential Analysis at Granscript-level resolution' Yi et al. 2017

# get all the transcript IDs
abundance <- read.table("data/fullTranscript/SRR1586618/abundance.tsv", header=T, sep="\t")
transcripts <- as.character(abundance$target_id)

transcript_split <- strsplit(transcripts, "_")
df_transcripts <- as.data.frame(do.call(rbind,transcript_split))
colnames(df_transcripts) <- c("gene", "transcript_id")

df_transcripts$target_id <- abundance$target_id

merge_results <- merge(sleuth_table, df_transcripts, by='target_id', all.x=TRUE)
head(merge_results)

source('analysis/scripts/aggregation.R')

MinMethod <- function(pvalues)
{
	if(length(pvalues) == 0)
	{
			return(NA)
	}
	pvalues <- pvalues[!is.na(pvalues)]
	if(length(pvalues) == 0)
	{
			return(NA)
	}
	n <- length(pvalues)
	m <- min(pvalues)
	1 - (1-m) ^ n
}

lancaster <- function(pvalues, weights)
{
	weights <- weights[!is.na(pvalues)]
	pvalues <- pvalues[!is.na(pvalues)]
	pvalues <- pvalues[weights>0]
	weights <- weights[weights>0]
	
	if(length(weights) != length(pvalues))
	{
		print('error, weights not equal to pvalues')
	}

	if(length(pvalues) == 0)
	{
		return(NA)
	}
	if(!any(weights))
	{
		return(NA)
	}
	if(length(pvalues) == 1)
	{
		return(pvalues)
	}
	t <- sapply(1:length(pvalues), function(i) lts(pvalues[i], weights[i]))
	t <- sum(t)
	p <- pchisq(t, sum(weights), lower.tail=FALSE) 
	return(p)
}

results <- merge_results %>% group_by(gene) %>% summarise(lan = lancaster(pval, exp(mean_obs)), weight = sum(exp(mean_obs), na.rm=TRUE), min = MinMethod(pval))
dim(results)

sig_results <- dplyr::filter(results, lan <= 0.05)
dim(sig_results)
head(sig_results)

plot_bootstrap(leaf_so, "GRMZM5G801958_T01", units="est_counts", color_by="cultivar")