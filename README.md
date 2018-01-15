# Highland Maize Genetic Incompatibility 

* Created: January 5, 2018
* Author: Caryn Johansen
* Contact: caryn.k.johansen@gmail.com

Last updated: January 15, 2018

Here, I will be going into more detail on the contents of the project directories in Projects/

**IMPORTANT NOTE:** While summarizing and cleaning my farm, I did move some files and directories. This may have changed how some of the files run - but I did not go in and change those scripts. If you are trying to run a script as is and it fails, the first issue checked should be the directory paths.

**MEA CULPA:** I wrote this all in vim very quickly. I did not spell check or edit.

## Info on Runcie Panel

https://github.com/deruncie/F1_selection/blob/master/Accession_selection_for_F1s.csv

## directory structure

## Highland_GeneticIncompatibility/

The parent proejct folder for the gene incompatibility project.

Highland_GeneticIncompatibility/
├── data
│   ├── processed
│   └── raw
├── LemmonAnalysis
│   ├── analysis
│   │   └── complete
│   ├── data
│   │   ├── processed
│   │   │   ├── v3
│   │   │   │   ├── bootstrap100
│   │   │   │   └── PRJNA262181_leaf
│   │   │   └── v4
│   │   │       ├── bootstrap100
│   │   │       ├── fullTranscript
│   │   │       └── PRJNA262181_leaf
│   │   └── raw
│   ├── logs
│   └── scripts
├── logs
├── scripts
└── XPCLR

Most of the current analysis in int he LemmonAnalysis directory, as that was the current focus of attention. As such, I will deal with that directory separately.

### data/

data/
├── processed
└── raw
    └── SRR1586618_pass_1.fastq.gz

The typical raw/ and processed/ directory structure. mostly Empty for now, except one test fastq in the raw/ directory.

### logs/

empty, but intended to be the output of slurm logs

### scripts

scripts/
└── fieldsummary.R

The only script for this is an R script meant to understand the number of plants that Dan Runcie had growing, and perhaps to eventually graph them. Data source in in the R script.

### XPCLR/

XPCLR/
├── MexLow_GuaHigh.allChr.wtclr.txt
├── MexLow_MexHigh.allChr.wtclr.txt
├── MexLow_SW_US.allChr.wtclr.txt
└── SA_Low_Andes.allChr.wtclr.txt

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5531767/

Updated XP-CLR values from above paper, supplied by Li Wang, for version 4 for the genome. Downloaded from iplant

## LemmonAnalisys/

There is a README.md file in LemmonAnalysis/ that going into detail for all the files, and provides data links, etc.

Here, I will provide an overview of the directories. This information is also in the LemmonAnalysis/README.md file.

Directory structure:

Highland_GeneticIncompatibility/LemmonAnalysis/
├── analysis
│   ├── GOanalysis
│   ├── images
│   ├── plotting_files
│   │   └── figure-html
│   └── scripts
├── data
│   ├── LemmonSupp
│   ├── processed
│   │   ├── v3
│   │   │   └── PRJNA262181_leaf
│   │   └── v4
│   │       └── PRJNA262181_leaf
│   └── raw
├── logs
├── paper_info
├── python_testing
├── run_info
└── scripts

### analysis/

The directory for the scripts analyzing the processed data, and for holding the outputs.


the .rds files are the saved outputs and steps of the R scripts.

* B73_TIL01_lfc.csv             - Lof fold change output between B73 and TIL01
* GOanalysis/                   - GO analysis R script
* images/                       - PDF and PNG files output from R scripts
* plotting.md                   - markdown output of plotting.Rmd
* plotting_files/               - image files for plotting.md
* scripts/                      - R scripts for analysis

#### analysis/scripts/

* aggregation.R                 - code from Yi et al. 2017 to aggregate the output p-value for kallisto and sleuth analysis
* basic_analysis.R              - some preliminary looks at kallisto abundance data
* plotting.Rmd                  - the main R script for plotting the output
* sleuthDE.R                    - test R code for using Sleuth

#### analysis/GOanalysis/

* maizeGOenrighment.R           - R code for conducting a GO enrichment test on lists of maize genes

#### analysis/images/

Image output from plotting.Rmd file

### data/

Highland_GeneticIncompatibility/LemmonAnalysis/data/
├── LemmonSupp
│   ├── SupplementalDataset1_v3.csv
│   └── SupplementalDataset1_v3.txt
├── processed
│   ├── v3
│   │   └── PRJNA262181_leaf
│   └── v4
│       └── PRJNA262181_leaf
└── raw

* LemmonSupp/                   - supplementary data set from Lemmon et al.
* processed/                    - kallisto output using v3 and v4 maize genomes
* processed/v3/PRJNA262181_leaf - kallisto output for the PRJNA262181 leaf samples for v3 maize genome
* processed/v4/PRJNA262181_leaf - kallisto output for the PRJNA262181 leaf samples for v4 maize genome
* raw/                          - fastq.gz files from SRA

### logs/

Slurm script error and standard outs

### paper_info/

Supplementary information from Lemmon et al. outlining the cross information.

### python_testing/

LemmonAnalysis/python_testing/
├── masterlist.txt
└── parse_kallisto_output.py

Python script to parse the kallisto output.

### run_info/

.info file as outputs of slurm scripts. Useful for debugging.


### scripts/

scripts/
├── kallisto_index.sh
├── kallisto.sh
├── kallisto_v3.sh
├── single_kallisto.sh
└── sra-dump.sh

* kallisto_index.sh                     - build the index for the reference genome using kallisto
* kallisto.sh                           - quantification of reads with kallisto using reference genome version 4
* kallisto_v3.sh                        - using maize reference genome version 3 for kallisto quantification
* sra-dump.sh                           - get data from SRA using samtools


