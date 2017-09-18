#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest
#SBATCH -J SRR1586620 
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

srr=SRR1586620

idx=$srr.idx
fasta=$srr\_pass.fasta.gz
fastq=$srr\_pass\_1.fastq.gz

## MAIN

echo $SLURM_JOB_ID $srr $fasta $fastq $idx >> $srr\_info.txt

#Build index

#idx_start=`date +%s`

#kallisto index -i data/raw/$idx data/raw/$fasta
#err=$?
#echo index error: $err

#idx_end=`date +%s`
#((idx_time=$end-$start))

#echo kallisto index runtime was $idx_time s

# Quantify

quant_start=`date +%s`

mkdir data/processed/$srr

kallisto quant -i data/raw/$idx -o data/processed/$srr -b 100 --single -l 101 -s 20 data/raw/$fastq
err=$?
echo kallisto error: $err

quant_end=`date +%s`
((quant_time=$quant_end - $quant_start))
echo kallisto quantification run time was $quant_time s

