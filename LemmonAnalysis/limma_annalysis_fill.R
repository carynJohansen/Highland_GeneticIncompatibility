##------------
# Libraries

library(limma)
library(tidyverse)
library(edgeR)
library(dplyr)

##------------
## wd

#setwd("/Users/caryn/Box Sync/Projects/Highland_GeneticIncompatibility/LemmonAnalysis")

#farm
setwd("/home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis")

##------------
# Functions

# function to read in a sample name and return the abundance file

get_abundance <- function(sample) {
	#sample is a character string
	#dir_path <- file.path(paste("data/PRJNA262181_leaf_v3/", sample, sep=""))
  dir_path_farm <- file.path(paste("data/processed/v3/PRJNA262181_leaf/", sample, sep=""))
	abundance_file <- paste(dir_path_farm, "abundance.tsv", sep="/")
	ab <- read.table(abundance_file, header=T, sep="\t")

	# need to separate the transcript ID and the gene name:
	transcripts <- as.character(ab$target_id)
	transcript_split <- strsplit(transcripts, "_")
	transcripts <- as.data.frame(do.call(rbind, transcript_split))
	colnames(transcripts) <- c("gene", "transcript_id")
	ab$gene <- transcripts$gene
	ab$transcript_id <- transcripts$transcript_id
	return(ab)
}


#get the tpm for all genes across transcripts
get_sum_tpm <- function(ab) {
	# where ab is the abundance data frame from get_abundance()
	# sum tpm over all transcripts of genes
	genes <- ab %>%
		group_by(gene) %>%
		summarise(est_tpm_sum=sum(tpm))
	genes <- data.frame(genes)
	colnames(genes) <- c("gene", "est_tpm_sum")

	return(genes)
}

# function for identifying the genotype, tissue, etc of sample name

get_sample_info <- function(sample) {
	leaf_sample_info <- leaf_list[[sample]]
	return(leaf_sample_info)
}

# function to generate the tpm matrix for edgeR DGE analysis
get_tpm_list <- function(v) {
  # to be used in an lapply function - v is a named vector of any length were the value named "sample" is the SRR number. No other values matter here.
  sample <- v["sample"]
  print(sample)
  ab <- get_abundance(sample)
  g_tpm <- get_sum_tpm(ab)
  rownames(g_tpm) <- g_tpm$gene
  g_tpm$gene <- NULL; colnames(g_tpm) <- sample
  #ab_list[[sample]] <- g_tpm
  return(g_tpm)
}

# function to get the DGE genes from edgeR

get_top <- function(cont_group, gene_tpm) {
  #where cont group should be a vector of three chacacters, gene_tpm is the full dataframe of all samples and their gene counts
  grp <- cont_group[1]
  print(grp)
  #get all sample names associated with the genotypes of the contrast group and subset the gene_tpm data.frame
  samples <- as.character(meta_leaf[meta_leaf$cultivar %in% cont_group,]$sample)
  print(samples)
  tpm_samples <- gene_tpm[, colnames(gene_tpm) %in% samples]
  
  # make a DGEList object
  dge <- DGEList(counts = tpm_samples, group = as.factor(samples))
  
  # calculate the form factors with calcNormFactors
  dge <- calcNormFactors(dge)
  
  # make a full contrast model
  sample_meta <- meta_leaf[meta_leaf$cultivar %in% cont_group,]
  sample_meta$cultivar <- factor(sample_meta$cultivar, levels=cont_group)
  sample_full <- model.matrix(~cultivar, data=sample_meta)
  
  # calculate the fit with voom, lmFit, and eBayes
  tpm_voom <- voom(dge, design = sample_full, plot=FALSE)
  fit <- lmFit(tpm_voom, design = sample_full)
  fit <- eBayes(fit)
  
  # calculate the maize parent contrast
  parent_maize <- cont_group[2]; parent_teo <- cont_group[3]
  prt_maize_crst <- paste("cultivar", parent_maize, sep="")
  prt_teo_crst <- paste("cultivar", parent_teo, sep="")
  
  fit_maize <- topTable(fit, n=dim(tpm_voom)[1], confint = TRUE, coef = prt_maize_crst)
  fit_maize$gene <- rownames(fit_maize)
  # calculate the teosinte parent contrast
  fit_teo <- topTable(fit, n=dim(tpm_voom)[1], confint = TRUE, coef = prt_teo_crst)
  fit_teo$gene <- rownames(fit_teo)
  
  # return some huge list
  res <- list(group = grp, tpm_samples = tpm_samples, samples = samples, contrasts = list(sample_full, prt_maize_crst, prt_teo_crst), 
              topMaize=fit_maize, topTeo=fit_teo)
}

# function to organize the results

get_maize_result <- function(res) {
  grp <- res[["group"]]
  maize_fit <- res[["topMaize"]]
  maize_fit$group <- grp
  return(maize_fit)
}
get_teo_result <- function(res) {
  grp <- res[["group"]]
  teo_fit <- res[["topTeo"]]
  teo_fit$group <- grp
  return(teo_fit)
}

# ploting tpm function
## NEEDS LONG DATA FRAME
plot_gene_tpm <- function(gene_to_plot) {
  #gene is a character string
  gene_tpm_long %>% filter(gene == gene_to_plot) %>% ggplot() + 
    geom_bar(aes(x=sample, y=tpm, fill=cultivar), stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(gene_to_plot)
}
# main function, to lapply over list of samples

##------------
# Data

# meta data linking sample names to gentoype and tissue
meta <- read.table("SraRunTable.txt", header=T, sep="\t")

# filter and organize the meta data
meta_leaf <- dplyr::filter(meta, tissue_s == "leaf")
meta_leaf <- dplyr::select(meta_leaf, sample=Run_s, cultivar=cultivar_s, genotype=genotype_s, tissue=tissue_s)
meta_leaf$cultivar <- as.character(meta_leaf$cultivar)
rm(meta)

# list of sample names, from directory or from meta_leaf

leaf_list <- list()
for (i in 1:dim(meta_leaf)[1]) {
	leaf_list[[i]] <- c(sample=as.character(meta_leaf[i,]$sample), genotype=as.character(meta_leaf[i,]$cultivar))
}
names(leaf_list) <- as.character(meta_leaf$sample)
leaf_list[[1]]
# contrasts
contrast_groups <- list("B73_TIL01"=c("B73_TIL01","B73", "TIL01"),
                        "B73_TIL03"=c("B73_TIL03","B73", "TIL03"),
                        "B73_TIL05"=c("B73_TIL05","B73", "TIL05"),
                        "B73_TIL09"=c("B73_TIL09","B73", "TIL09"),
                        "B73_TIL11"=c("B73_TIL11","B73", "TIL11"),
                        "B73_TIL14"=c("B73_TIL14","B73", "TIL14"),
                        "B73_TIL25"=c("B73_TIL25","B73", "TIL25"),
                        "CML103_TIL03"=c("CML103_TIL03","CML103", "TIL03"),
                        "CML103_TIL11"=c("CML103_TIL11","CML103", "TIL11"),
                        "CML103_TIL14"=c("CML103_TIL14","CML103", "TIL14"),
                        "CML103_TIL25"=c("CML103_TIL25","CML103", "TIL25"),
                        "Ki3_TIL03"=c("Ki3_TIL03","Ki3", "TIL03"),
                        "Ki3_TIL09"=c("Ki3_TIL09","Ki3", "TIL09"),
                        "Ki3_TIL11"=c("Ki3_TIL11","Ki3", "TIL11"),
                        "Ki3_TIL14"=c("Ki3_TIL14","Ki3", "TIL14"),
                        "Mo17_TIL01"=c("Mo17_TIL01","Mo17", "TIL01"),
                        "Mo17_TIL09"=c("Mo17_TIL09","Mo17", "TIL09"),
                        "Mo17_TIL14"=c("Mo17_TIL14","Mo17", "TIL14"),
                        "Mo17_TIL25"=c("Mo17_TIL25","Mo17", "TIL25"),
                        "Oh43_TIL01"=c("Oh43_TIL01","Oh43", "TIL01"),
                        "Oh43_TIL03"=c("Oh43_TIL03","Oh43", "TIL03"),
                        "Oh43_TIL09"=c("Oh43_TIL09","Oh43", "TIL09"),
                        "Oh43_TIL10"=c("Oh43_TIL10","Oh43", "TIL10"),
                        "Oh43_TIL11"=c("Oh43_TIL11","Oh43", "TIL11"),
                        "Oh43_TIL15"=c("Oh43_TIL15","Oh43", "TIL15"),
                        "Oh43_TIL25"=c("Oh43_TIL25","Oh43", "TIL25"),
                        "W22_TIL01"=c("W22_TIL01","W22", "TIL01"),
                        "W22_TIL03"=c("W22_TIL03","W22", "TIL03"),
                        "W22_TIL11"=c("W22_TIL11","W22", "TIL11"),
                        "W22_TIL14"=c("W22_TIL14","W22", "TIL14"),
                        "W22_TIL25"=c("W22_TIL25","W22", "TIL25"))#,

parent_contrasts <- list(
  "B73_CM103"=c("B73","CM103"),
  "B")


##------------
# Main

# for each sample, get the tpm for each gene, put into a matrix where rows are gene, columns are sample
ab_list <- lapply(leaf_list, get_tpm_list)
gene_tpm <- do.call(cbind, ab_list)
saveRDS(gene_tpm, file="gene_tpm.rds")
# topTables lists for each contrast group
result_list <- lapply(contrast_groups, get_top, gene_tpm)
saveRDS(result_list, file="full_result_list.rds")

maize_results_list <- lapply(result_list, get_maize_result)
maize_results <- do.call(rbind, maize_results_list)

teo_results_list <- lapply(result_list, get_teo_result)
teo_results <- do.call(rbind, teo_results_list)

saveRDS(maize_results, file="maize_topTable_results.rds")
saveRDS(teo_results, file="teo_topTable_results.rds")

maize_results <- readRDS("maize_topTable_results.rds")
teo_results <- readRDS("teo_topTable_results.rds")
# plot all
res_all <- merge(maize_results, teo_results, by=c("gene","group"))
saveRDS(res_all, "merged_results.rds")

res_all$group <- as.factor(res_all$group)


png(file="all_results.png")
ggplot(res_all) + geom_point(aes(x=logFC.x, y=logFC.y, color=group), alpha=0.5)
dev.off()

# significant genes in both
res_sig <- res_all %>% filter(adj.P.Val.x < 0.001 | adj.P.Val.y < 0.001)

res_sig <- res_sig %>% mutate(de = ifelse(adj.P.Val.x < 0.001 & adj.P.Val.y > 0.001, "maize",
  ifelse(adj.P.Val.x > 0.001 & adj.P.Val.y < 0.001, "teosinte",
    "both")))

dim(res_sig)
table(res_sig$group)
table(res_sig$de)

saveRDS(res_sig, file="bygroup_sig.rds")

pdf(file="sig_genes.pdf")
ggplot(res_sig) + geom_point(aes(x=logFC.x, y=logFC.y, color=group), alpha=0.5)
dev.off()

png(file="sig_genes.png")
ggplot(res_sig) + geom_point(aes(x=logFC.x, y=logFC.y, color=group), alpha=0.5)
dev.off()

# get the unique genes that are significant for both
res_sig_both  <- res_sig %>% filter(de == "both")
N_unique_genes <- n_distinct(res_sig_both$gene)


png(file="sig_genes_both.png")
ggplot(res_sig_both) + geom_point(aes(x=logFC.x, y=logFC.y, color=group), alpha=0.5)
dev.off()
pdf(file="sig_genes_both.pdf")
ggplot(res_sig_both) + geom_point(aes(x=logFC.x, y=logFC.y, color=group), alpha=0.5)
dev.off()
# genes under/over dom common/shared in maize

# genes under/over dom common/shared in teosinte

# remove those genes, find genes just under/over dom in F1s

# plot those genes

##-------------------
## Genotype Contrasts

# I'm conserned about the multiple testing here.
# so an alterative idea, rather than testing by cultivar, use genotype ("maize", "teosinte", "F1 cross")

genotype_full <- model.matrix(~genotype, data=meta_leaf)
genotype_dge <- DGEList(counts=gene_tpm, group=meta_leaf$genotype)
genotype_dge <- calcNormFactors(genotype_dge)

genotype_voom <- voom(genotype_dge, 
  design = genotype_full, plot=FALSE)
genotype_fit <- lmFit(genotype_voom, design = genotype_full)
genotype_fit <- eBayes(genotype_fit)

genotype_top <- topTable(genotype_fit, n=Inf, confint = TRUE)
dim(genotype_top)
saveRDS(genotype_top, file="byGenotype_top.rds")

png(file="genotype_all.png")
ggplot(genotype_top) + geom_point(aes(x=genotypemaize, y=genotypeteosinte), alpha=0.5)
dev.off()

genotype_sig <- genotype_top %>% filter(adj.P.Val < 0.001)
dim(genotype_sig)

png(file="genotype_sig.png")
ggplot(genotype_sig) + geom_point(aes(x=genotypemaize, y=genotypeteosinte), alpha=0.5)
dev.off()

# Calculate for maize and teosinte specifically

maize_geno_top <- topTable(genotype_fit, n=Inf, confint=TRUE, coef = "genotypemaize")
teo_geno_top <- topTable(genotype_fit, n=Inf, confint=TRUE, coef = "genotypeteosinte")

# get sig genes for each comparison
maize_geno_sig <- maize_geno_top %>% filter(adj.P.Val < 0.001)
teo_geno_sig <- teo_geno_top %>% filter(adj.P.Val < 0.001)
dim(teo_geno_sig)

# merge into one data set

maize_geno_top$gene <- rownames(maize_geno_top)
teo_geno_top$gene <- rownames(teo_geno_top)
geno_all <- merge(maize_geno_top, teo_geno_top, by="gene")
head(geno_all)
geno_all_sig <- geno_all %>% filter(adj.P.Val.x < 0.001 | adj.P.Val.y < 0.001)
dim(geno_all_sig)
# make a column for which genotype was significant, maize, teosinte, or both

geno_all_sig <- geno_all_sig %>% mutate(de = ifelse(adj.P.Val.x < 0.001 & adj.P.Val.y > 0.001, "maize",
  ifelse(adj.P.Val.x > 0.001 & adj.P.Val.y < 0.001, "teosinte",
    "both")))
geno_all_sig$de <- as.factor(as.character(geno_all_sig$de))
saveRDS(geno_all_sig, file="bygenotype_sig.rds")

teo_sig <- geno_all_sig %>% filter(adj.P.Val.y < 0.001 & adj.P.Val.x > 0.001)
dim(teo_sig)
pdf(file="merged_geno_sig.pdf")
ggplot(geno_all_sig) + geom_point(aes(x=logFC.x, y=logFC.y, color=de), alpha = 0.5)
dev.off()

##----------------
## comaprison to mid-parent

# calculate value you might expect based on mid-parent
# test for difference from mid-parent

##----------------
## Ploting gene tpm

library(reshape2)
gene_tpm_long <- gene_tpm
gene_tpm_long$gene <- rownames(gene_tpm)
gene_tpm_long <- reshape2::melt(gene_tpm_long)
colnames(gene_tpm_long) <- c("gene", "sample", "tpm")
gene_tpm_long <- gene_tpm_long %>% left_join(meta_leaf, by = "sample")