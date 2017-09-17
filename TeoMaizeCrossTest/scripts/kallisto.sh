#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest
#SBATCH -J kallisto
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/kallisto_%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/kallisto_%j.out
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

outdir="/home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/data/processed"

idx=$srr.idx
fasta=$srr\_pass.fasta.gz
fastq=$srr\_pass\_1.fastq.gz

## MAIN

echo $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID $srr $fasta $fastq $idx >> kallisto_info.txt

#Build index

idx_start=`date +%s`

kallisto index -i data/raw/$idx data/raw/$fasta
err=$?
echo index error: $err

idx_end=`date +%s`
((idx_time=$end-$start))

echo kallisto index runtime was $idx_time s

# Quantify

quant_start=`date +%s`

mkdir data/processed/$srr

kallisto quant -i data/raw/$idx -o data/processed/$srr -b 100 --single -l 101 -s 20 data/raw/$fastq
err=$?
echo kallisto error: $err

quant_end=`date +%s`
((quant_time=$quant_end - $quant_start))
echo kallisto quantification run time was $quant_time s

