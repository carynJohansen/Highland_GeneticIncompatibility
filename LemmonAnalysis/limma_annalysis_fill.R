##------------
# Libraries

library(limma)
library(tidyverse)
library(edgeR)

##------------
## wd

setwd("/Users/caryn/Box Sync/Projects/Highland_GeneticIncompatibility/LemmonAnalysis")

#farm
#setwd("/home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis")

##------------
# Functions

# function to read in a sample name and return the abundance file

get_abundance <- function(sample) {
	sample <- "SRR1586766"
	#sample is a character string
	dir_path <- file.path(paste("data/PRJNA262181_leaf_v3/", sample, sep=""))
	abundance_file <- paste(dir_path, "abundance.tsv", sep="/")
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
  ab <- get_abundance(sample)
  g_tpm <- get_sum_tpm(ab)
  rownames(g_tpm) <- g_tpm$gene
  g_tpm$gene <- NULL; colnames(g_tpm) <- sample
  #ab_list[[sample]] <- g_tpm
  return(g_tpm)
}

# function to get the DGE genes from edgeR

dge_list <- function(cont_group, gene_tpm) {
  #where cont group should be a vector of three chacacters, gene_tpm is the full dataframe of all samples and their gene counts
  grp <- cont_group[3]
  #get all sample names associated with the genotypes of the contrast group and subset the gene_tpm data.frame
  samples <- as.character(meta_leaf[meta_leaf$cultivar %in% cont_group,]$sample)
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
  parent_maize <- cont_group[1]; parent_teo <- cont_group[2]
  prt_maize_crst <- paste("cultivar", parent_maize, sep="")
  prt_teo_crst <- paste("cultivar", parent_teo, sep="")
  
  fit_maize <- topTable(fit, n=nrow(tpm_voom), confint = TRUE, coef = prt_maize_crst)
  fit_maize$gene <- rownames(fit_maize)
  # calculate the teosinte parent contrast
  fit_teo <- topTable(fit, n=nrow(tpm_voom), confint = TRUE, coef = prt_teo_crst)
  fit_teo$gene <- rownames(fit_teo)
  
  # return some huge list
  res <- list(group = grp, tpm_samples = tpm_samples, samples = samples, contrasts = list(sample_full, prt_maize_crst, prt_teo_crst), 
              topMaize=fit_maize, topTeo=fit_teo)
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
meta_leaf <- meta %>% filter(tissue_s == "leaf") %>% select(sample=Run_s, cultivar=cultivar_s, genotype=genotype_s, tissue=tissue_s)
meta_leaf$cultivar <- as.character(meta_leaf$cultivar)
rm(meta)

# list of sample names, from directory or from meta_leaf

leaf_list <- list()
for (i in 1:nrow(meta_leaf)) {
	leaf_list[[i]] <- c(sample=as.character(meta_leaf[i,]$sample), genotype=as.character(meta_leaf[i,]$cultivar))
}
names(leaf_list) <- as.character(meta_leaf$sample)

# contrasts
contrast_groups <- list("B73_TIL01"=c("B73", "TIL01", "B73_TIL01"),
                        "B73_TIL03"=c("B73", "TIL03", "B73_TIL03"),
                        "W22_TIL01"=c("W22", "TIL01", "B73_TIL03"))#,



##------------
# Main

# for each sample, get the tpm for each gene, put into a matrix where rows are gene, columns are sample
ab_list <- lapply(leaf_list, get_tpm_list)
gene_tpm <- do.call(cbind, ab_list)

# topTables lists for each contrast group
result_list <- lapply(contrast_groups, dge_list, gene_tpm)


# plot all

# genes under/over dom common/shared in maize

# genes under/over dom common/shared in teosinte

# remove those genes, find genes just under/over dom in F1s

# plot those genes