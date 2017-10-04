#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis
#SBATCH -J kallisto
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/T01_kallisto_%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/T01_kallisto_%j.out
#SBATCH --time=9:00:00
#SBATCH --mem=50000
#SBATCH --array=1-13

set -u
set -e

## MODULES

module load kallisto

#print version
kallisto version

## DATA

srrFile="srr_numbers.txt"

srr=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR == var {print}' $srrFile)


idx="/home/caryn89/genomes/maize_cdna_v4/Zea_mays.AGPv4.cdna.idx"

fastq=$srr\_pass\_1.fastq.gz

## MAIN

# Quantify

quant_start=`date +%s`


mkdir data/processed/v4/$srr

kallisto quant -i $idx -o data/processed/v4/$srr -b 50 --single -l 101 -s 20 data/raw/$fastq

err=$?
echo kallisto error: $err

quant_end=`date +%s`
((quant_time=$quant_end - $quant_start))
echo kallisto quantification run time was $quant_time s

echo $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID $quant_time $srr $fastq $idx >> b100_kallisto_info.txt

