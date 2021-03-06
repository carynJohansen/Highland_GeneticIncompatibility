#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/ 
#SBATCH -o /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/sra-leaves-dump-%j.out
#SBATCH -e /home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/logs/sra-leaves-dump-%j.out
#SBATCH -J fasta
#SBATCH --time=5:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cjohansen@ucdavis.edu
#SBATCH --mem=4000
#SBATCH --array=4,5,76,77,78,79,80,82,83

set -u
set -e

## MODULES

module load sratoolkit/2.8.2

## DATA

srrFile="leaf_srr_numbers.txt"

srr=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR == var {print}' $srrFile)
echo $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID $srr >> leaves_sra_info.txt

#for if
fastq=$srr\_pass\_1.fastq.gz

outdir="/home/caryn89/Projects/Highland_GeneticIncompatibility/LemmonAnalysis/data/raw"

#fastq-dump --outdir data/raw --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip $srr

#check to see if the file has already been created

if [ ! -e $outdir/$fastq ]; then
	#fastq-dump --outdir data/raw --gzip --skip-technical --readids --read-filter pass --dumpbase --clip -fasta $srr
	fastq-dump --outdir data/raw --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip $srr
	exitError=$?
	echo exit error $exitError
else
	echo $fastq "already exists"
	echo "located at: " $outdir
fi


