#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis
#SBATCH -J v3index
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/v3index_%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/v3index_%j.out
#SBATCH --time=5:00:00
#SBATCH --mem=50000

set -u
set -e

## MODULES

module load kallisto

#print version
kallisto version

## DATA

cdna="/home/caryn89/genomes/maize_cdna_v3/Zea_mays.AGPv3.cdna.all.fa.gz"

idx="/home/caryn89/genomes/maize_cdna_v3/Zea_mays.AGPv3.cdna.all.fa.idx"


## MAIN

#Build index

idx_start=`date +%s`

kallisto index -i $idx $cdna
err=$?
echo index error: $err

idx_end=`date +%s`
((idx_time=$idx_end-$idx_start))
echo $idx_time

