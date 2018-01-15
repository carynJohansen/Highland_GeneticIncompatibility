# Testing the Gene Incompatibility hypothesis/analysis pipeline on Maize:Teosinte data

The purpose of this example project is to use RNA-Seq data from the paper Lemmon et al. 2014.

Paper found here: http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004745

Lemmon et al. made crosses between teosinte inbred lines and maize inbred lines.

Download data from SRA:

SRA project: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA262181

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61810

meta data to be found at: https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_84857220_130.14.18.48_5555_1508787163_1278959485_0MetA0_S_HStore&query_key=2

# Analysis

For differential gene expression analysis, see limmaDE.Rmd

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
* sra-dump.sh				- get data from SRA using samtools

