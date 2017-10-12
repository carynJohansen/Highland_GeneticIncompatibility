#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis
#SBATCH -J leavesv3
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/v3leaf_kallisto_%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/v3leaf_kallisto_%j.out
#SBATCH --time=9:00:00
#SBATCH --mem=50000
#SBATCH --array=1-151

set -u
set -e

## MODULES

module load kallisto

#print version
kallisto version

## DATA

srrFile="leaf_srr_numbers.txt"

srr=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR == var {print}' $srrFile)

echo $srr

idx="/home/caryn89/genomes/maize_cdna_v3/Zea_mays.AGPv3.cdna.all.fa.idx"

fastq=$srr\_pass\_1.fastq.gz

echo $fastq
echo $idx

## MAIN

# Quantify

quant_start=`date +%s`


mkdir data/processed/v3/PRJNA262181_leaf/$srr

kallisto quant -i $idx -o data/processed/v3/PRJNA262181_leaf/$srr -b 100 --single -l 101 -s 20 data/raw/$fastq

err=$?
echo kallisto error: $err

quant_end=`date +%s`
((quant_time=$quant_end - $quant_start))
echo kallisto quantification run time was $quant_time s

echo $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID $quant_time $srr >> PRJNA262181_leaf_kallisto_info.txt

