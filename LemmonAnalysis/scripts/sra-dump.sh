#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/ 
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/sra-dump-%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/logs/sra-dump-%j.out
#SBATCH -J fasta
#SBATCH --time=5:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cjohansen@ucdavis.edu
#SBATCH --mem=4000
#SBATCH --array=1-13

set -u
set -e

## MODULES

module load sratoolkit/2.8.2

## DATA

srrFile="srr_numbers.txt"

srr=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR == var {print}' $srrFile)
echo $SLURM_ARRAY_TASK_ID $srr >> sra_info.txt

outdir="/home/caryn89/Projects/Highland_GeneticIncompatibility/TeoMaizeCrossTest/data/raw"

#fastq-dump --outdir data/raw --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip $srr

fastq-dump --outdir data/raw --gzip --skip-technical --readids --read-filter pass --dumpbase --clip -fasta $srr

exitError=$?
echo exit error $exitError


