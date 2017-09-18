#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest
#SBATCH -J index
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/index_%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/index_%j.out
#SBATCH --time=5:00:00
#SBATCH --mem=50000

set -u
set -e

## MODULES

module load kallisto

#print version
kallisto version

## DATA

cdna="/home/caryn89/genomes/maize_cdna_v4/Zea_mays.AGPv4.cdna.all.fa.gz"

idx="/home/caryn89/genomes/maize_cdna_v4/Zea_mays.AGPv4.cdna.all.idx"


## MAIN

#Build index

idx_start=`date +%s`

kallisto index -i $idx $cdna
err=$?
echo index error: $err

idx_end=`date +%s`
((idx_time=$idx_end-$idx_start))
echo $index_time

