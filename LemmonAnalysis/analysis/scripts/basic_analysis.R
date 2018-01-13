# Clean the headers from the fasta to include only sequences with transcript id 01 or 001


setwd("/Users/caryn/Box Sync/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest")

library(tidyverse)
library(ggplot2)

#SRR1586618


abundance <- read.table("data/processed/SRR1586619/abundance.tsv", header=T, sep="\t")

#Filter data for only single transcripts

targets <- as.character(abundance$target_id)

tar_split <- strsplit(targets, "_")

tars <- do.call(rbind.data.frame, tar_split)
colnames(tars) <- c("geneID", "transcriptID")

abundance <- cbind(abundance, tars)

# Check it out

abundance %>%
  summarise(min_est_counts=min(abundance$est_counts), 
            max_est_counts=max(abundance$est_counts), 
            sd_est_counts=sd(abundance$est_counts),
            min_gene_length=min(abundance$length),
            max_gene_length=max(abundance$length))

max <- which(abundance$est_counts == max(abundance$est_counts))
abundance[max,]

abundance %>%
  ggplot() + geom_density(aes(x=log(est_counts)))

abundance %>%
  ggplot() + geom_histogram(aes(x=est_counts),bins = 100)

abundance %>% 
  ggplot() + geom_point(aes(x=tpm, y=est_counts), alpha = 1/10)

abundance %>%
  ggplot() + geom_point(aes(x=length, y=est_counts), alpha = 1/8)

abundance %>% 
  ggplot() + geom_bar(aes(x=transcriptID, y=est_counts), stat = "identity")

abundance %>% 
  ggplot() + geom_bar(aes(x=transcriptID, y=tpm), stat = "identity")


