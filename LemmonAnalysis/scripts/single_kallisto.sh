#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest
#SBATCH -J SRR1586766 
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/single_kal_%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/single_kal_%j.out
#SBATCH --time=9:00:00
#SBATCH --mem=50000

set -u
set -e

## MODULES

module load kallisto

#print version
kallisto version

## DATA

srr=SRR1586766

idx="/home/caryn89/genomes/maize_cdna_v4/Zea_mays.AGPv4.cdna.all.idx"
fastq=$srr\_pass\_1.fastq.gz

## MAIN

echo $SLURM_JOB_ID $srr $fastq $idx >> $srr\_info.txt

# Quantify

quant_start=`date +%s`

mkdir data/processed/$srr

#not bootstrapping for test runs
kallisto quant -i $idx -o data/processed/$srr -b 100 --single -l 101 -s 20 data/raw/$fastq
err=$?
echo kallisto error: $err

quant_end=`date +%s`
((quant_time=$quant_end - $quant_start))
echo kallisto quantification run time was $quant_time s

