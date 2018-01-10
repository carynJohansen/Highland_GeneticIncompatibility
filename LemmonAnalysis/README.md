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

Directory overview:

* all_srr_numbers.txt		==> list of SRR numbers from Lemmon et al. \n separated
* analysis/			==> parent directory of outputs/
* b100_kallisto_info.txt	==> bash output from script: scripts/kallk

Directory details:

LemmonAnalysis/
├── all_srr_numbers.txt
├── analysis
│   ├── bygenotype_sig.rds
│   ├── bygroup_sig.rds
│   ├── complete
│   │   ├── all_results.pdf
│   │   ├── all_results.png
│   │   ├── full_result_list.rds
│   │   ├── gene_tpm.rds
│   │   ├── genotype_all.png
│   │   ├── genotype_sig.png
│   │   ├── maize_topTable_results.rds
│   │   ├── merged_geno_sig.pdf
│   │   ├── merged_results.rds
│   │   ├── sig_genes.pdf
│   │   ├── sig_genes.png
│   │   └── teo_topTable_results.rds
│   ├── leaf_sleuth_results.rds
│   └── leaf_table.rds
├── b100_kallisto_info.txt
├── data
│   ├── processed
│   │   ├── v3
│   │   │   ├── bootstrap100
│   │   │   │   ├── SRR1586618
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586619
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586620
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586621
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586736
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586747
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586766
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586767
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586768
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586888
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586889
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   ├── SRR1586898
│   │   │   │   │   ├── abundance.h5
│   │   │   │   │   ├── abundance.tsv
│   │   │   │   │   └── run_info.json
│   │   │   │   └── SRR1586899
│   │   │   │       ├── abundance.h5
│   │   │   │       ├── abundance.tsv
│   │   │   │       └── run_info.json
│   │   │   └── PRJNA262181_leaf
│   │   │       ├── SRR1586766
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586767
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586768
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586769
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586770
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586771
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586772
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586773
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586774
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586775
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586776
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586777
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586778
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586779
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586780
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586781
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586782
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586783
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586784
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586785
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586786
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586787
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586788
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586789
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586790
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586791
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586792
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586793
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586794
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586795
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586796
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586797
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586798
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586799
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586800
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586801
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586802
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586803
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586804
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586805
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586806
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586807
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586808
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586809
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586810
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586811
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586812
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586813
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586814
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586815
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586816
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586817
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586818
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586819
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586820
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586821
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586822
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586823
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586824
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586825
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586826
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586827
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586828
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586829
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586830
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586831
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586832
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586833
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586834
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586835
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586836
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586837
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586838
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586839
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586840
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586841
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586842
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586843
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586844
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586845
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586846
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586847
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586848
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586849
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586850
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586851
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586852
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586853
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586854
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586855
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586856
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586857
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586858
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586859
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586860
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586861
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586862
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586863
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586864
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586865
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586866
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586867
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586868
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586869
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586870
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586871
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586872
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586873
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586874
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586875
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586876
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586877
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586878
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586879
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586880
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586881
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586882
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586883
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586884
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586885
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586886
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586887
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586888
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586889
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586890
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586891
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586892
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586893
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586894
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586895
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586896
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586897
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586898
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586899
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586900
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586901
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586902
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586903
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586904
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586905
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586906
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586907
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586908
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586909
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586910
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586911
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586912
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586913
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586914
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       ├── SRR1586915
│   │   │       │   ├── abundance.h5
│   │   │       │   ├── abundance.tsv
│   │   │       │   └── run_info.json
│   │   │       └── SRR1586916
│   │   │           ├── abundance.h5
│   │   │           ├── abundance.tsv
│   │   │           └── run_info.json
│   │   └── v4
│   │       ├── bootstrap100
│   │       │   ├── SRR1586618
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586619
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586620
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586621
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586736
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586747
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586766
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586767
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586768
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586888
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586889
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586898
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   └── SRR1586899
│   │       │       ├── abundance.h5
│   │       │       ├── abundance.tsv
│   │       │       └── run_info.json
│   │       ├── fullTranscript
│   │       │   ├── SRR1586618
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586619
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586620
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586621
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586736
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586747
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586766
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586767
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586768
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586888
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586889
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   ├── SRR1586898
│   │       │   │   ├── abundance.h5
│   │       │   │   ├── abundance.tsv
│   │       │   │   └── run_info.json
│   │       │   └── SRR1586899
│   │       │       ├── abundance.h5
│   │       │       ├── abundance.tsv
│   │       │       └── run_info.json
│   │       └── PRJNA262181_leaf
│   │           ├── SRR1586766
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586767
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586768
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586769
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586770
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586771
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586772
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586773
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586774
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586775
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586776
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586777
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586778
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586779
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586780
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586781
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586782
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586783
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586784
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586785
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586786
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586787
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586788
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586789
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586790
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586791
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586792
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586793
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586794
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586795
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586796
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586797
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586798
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586799
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586800
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586801
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586802
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586803
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586804
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586805
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586806
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586807
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586808
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586809
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586810
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586811
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586812
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586813
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586814
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586815
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586816
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586817
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586818
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586819
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586820
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586821
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586822
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586823
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586824
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586825
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586826
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586827
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586828
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586829
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586830
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586831
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586832
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586833
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586834
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586835
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586836
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586837
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586838
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586839
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586840
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586841
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586842
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586843
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586844
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586845
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586846
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586847
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586848
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586849
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586850
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586851
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586852
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586853
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586854
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586855
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586856
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586857
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586858
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586859
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586860
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586861
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586862
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586863
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586864
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586865
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586866
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586867
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586868
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586869
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586870
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586871
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586872
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586873
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586874
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586875
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586876
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586877
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586878
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586879
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586880
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586881
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586882
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586883
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586884
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586885
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586886
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586887
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586888
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586889
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586890
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586891
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586892
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586893
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586894
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586895
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586896
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586897
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586898
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586899
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586900
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586901
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586902
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586903
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586904
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586905
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586906
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586907
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586908
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586909
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586910
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586911
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586912
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586913
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586914
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           ├── SRR1586915
│   │           │   ├── abundance.h5
│   │           │   ├── abundance.tsv
│   │           │   └── run_info.json
│   │           └── SRR1586916
│   │               ├── abundance.h5
│   │               ├── abundance.tsv
│   │               └── run_info.json
│   └── raw
│       ├── SRR1586618_pass_1.fastq.gz
│       ├── SRR1586619_pass_1.fastq.gz
│       ├── SRR1586620_pass_1.fastq.gz
│       ├── SRR1586621_pass_1.fastq.gz
│       ├── SRR1586736_pass_1.fastq.gz
│       ├── SRR1586747_pass_1.fastq.gz
│       ├── SRR1586766_pass_1.fastq.gz
│       ├── SRR1586767_pass_1.fastq.gz
│       ├── SRR1586768_pass_1.fastq.gz
│       ├── SRR1586769_pass_1.fastq.gz
│       ├── SRR1586770_pass_1.fastq.gz
│       ├── SRR1586771_pass_1.fastq.gz
│       ├── SRR1586772_pass_1.fastq.gz
│       ├── SRR1586773_pass_1.fastq.gz
│       ├── SRR1586774_pass_1.fastq.gz
│       ├── SRR1586775_pass_1.fastq.gz
│       ├── SRR1586776_pass_1.fastq.gz
│       ├── SRR1586777_pass_1.fastq.gz
│       ├── SRR1586778_pass_1.fastq.gz
│       ├── SRR1586779_pass_1.fastq.gz
│       ├── SRR1586780_pass_1.fastq.gz
│       ├── SRR1586781_pass_1.fastq.gz
│       ├── SRR1586782_pass_1.fastq.gz
│       ├── SRR1586783_pass_1.fastq.gz
│       ├── SRR1586784_pass_1.fastq.gz
│       ├── SRR1586785_pass_1.fastq.gz
│       ├── SRR1586786_pass_1.fastq.gz
│       ├── SRR1586787_pass_1.fastq.gz
│       ├── SRR1586788_pass_1.fastq.gz
│       ├── SRR1586789_pass_1.fastq.gz
│       ├── SRR1586790_pass_1.fastq.gz
│       ├── SRR1586791_pass_1.fastq.gz
│       ├── SRR1586792_pass_1.fastq.gz
│       ├── SRR1586793_pass_1.fastq.gz
│       ├── SRR1586794_pass_1.fastq.gz
│       ├── SRR1586795_pass_1.fastq.gz
│       ├── SRR1586796_pass_1.fastq.gz
│       ├── SRR1586797_pass_1.fastq.gz
│       ├── SRR1586798_pass_1.fastq.gz
│       ├── SRR1586799_pass_1.fastq.gz
│       ├── SRR1586800_pass_1.fastq.gz
│       ├── SRR1586801_pass_1.fastq.gz
│       ├── SRR1586802_pass_1.fastq.gz
│       ├── SRR1586803_pass_1.fastq.gz
│       ├── SRR1586804_pass_1.fastq.gz
│       ├── SRR1586805_pass_1.fastq.gz
│       ├── SRR1586806_pass_1.fastq.gz
│       ├── SRR1586807_pass_1.fastq.gz
│       ├── SRR1586808_pass_1.fastq.gz
│       ├── SRR1586809_pass_1.fastq.gz
│       ├── SRR1586810_pass_1.fastq.gz
│       ├── SRR1586811_pass_1.fastq.gz
│       ├── SRR1586812_pass_1.fastq.gz
│       ├── SRR1586813_pass_1.fastq.gz
│       ├── SRR1586814_pass_1.fastq.gz
│       ├── SRR1586815_pass_1.fastq.gz
│       ├── SRR1586816_pass_1.fastq.gz
│       ├── SRR1586817_pass_1.fastq.gz
│       ├── SRR1586818_pass_1.fastq.gz
│       ├── SRR1586819_pass_1.fastq.gz
│       ├── SRR1586820_pass_1.fastq.gz
│       ├── SRR1586821_pass_1.fastq.gz
│       ├── SRR1586822_pass_1.fastq.gz
│       ├── SRR1586823_pass_1.fastq.gz
│       ├── SRR1586824_pass_1.fastq.gz
│       ├── SRR1586825_pass_1.fastq.gz
│       ├── SRR1586826_pass_1.fastq.gz
│       ├── SRR1586827_pass_1.fastq.gz
│       ├── SRR1586828_pass_1.fastq.gz
│       ├── SRR1586829_pass_1.fastq.gz
│       ├── SRR1586830_pass_1.fastq.gz
│       ├── SRR1586831_pass_1.fastq.gz
│       ├── SRR1586832_pass_1.fastq.gz
│       ├── SRR1586833_pass_1.fastq.gz
│       ├── SRR1586834_pass_1.fastq.gz
│       ├── SRR1586835_pass_1.fastq.gz
│       ├── SRR1586836_pass_1.fastq.gz
│       ├── SRR1586837_pass_1.fastq.gz
│       ├── SRR1586838_pass_1.fastq.gz
│       ├── SRR1586839_pass_1.fastq.gz
│       ├── SRR1586840_pass_1.fastq.gz
│       ├── SRR1586841_pass_1.fastq.gz
│       ├── SRR1586842_pass_1.fastq.gz
│       ├── SRR1586843_pass_1.fastq.gz
│       ├── SRR1586844_pass_1.fastq.gz
│       ├── SRR1586845_pass_1.fastq.gz
│       ├── SRR1586846_pass_1.fastq.gz
│       ├── SRR1586847_pass_1.fastq.gz
│       ├── SRR1586848_pass_1.fastq.gz
│       ├── SRR1586849_pass_1.fastq.gz
│       ├── SRR1586850_pass_1.fastq.gz
│       ├── SRR1586851_pass_1.fastq.gz
│       ├── SRR1586852_pass_1.fastq.gz
│       ├── SRR1586853_pass_1.fastq.gz
│       ├── SRR1586854_pass_1.fastq.gz
│       ├── SRR1586855_pass_1.fastq.gz
│       ├── SRR1586856_pass_1.fastq.gz
│       ├── SRR1586857_pass_1.fastq.gz
│       ├── SRR1586858_pass_1.fastq.gz
│       ├── SRR1586859_pass_1.fastq.gz
│       ├── SRR1586860_pass_1.fastq.gz
│       ├── SRR1586861_pass_1.fastq.gz
│       ├── SRR1586862_pass_1.fastq.gz
│       ├── SRR1586863_pass_1.fastq.gz
│       ├── SRR1586864_pass_1.fastq.gz
│       ├── SRR1586865_pass_1.fastq.gz
│       ├── SRR1586866_pass_1.fastq.gz
│       ├── SRR1586867_pass_1.fastq.gz
│       ├── SRR1586868_pass_1.fastq.gz
│       ├── SRR1586869_pass_1.fastq.gz
│       ├── SRR1586870_pass_1.fastq.gz
│       ├── SRR1586871_pass_1.fastq.gz
│       ├── SRR1586872_pass_1.fastq.gz
│       ├── SRR1586873_pass_1.fastq.gz
│       ├── SRR1586874_pass_1.fastq.gz
│       ├── SRR1586875_pass_1.fastq.gz
│       ├── SRR1586876_pass_1.fastq.gz
│       ├── SRR1586877_pass_1.fastq.gz
│       ├── SRR1586878_pass_1.fastq.gz
│       ├── SRR1586879_pass_1.fastq.gz
│       ├── SRR1586880_pass_1.fastq.gz
│       ├── SRR1586881_pass_1.fastq.gz
│       ├── SRR1586882_pass_1.fastq.gz
│       ├── SRR1586883_pass_1.fastq.gz
│       ├── SRR1586884_pass_1.fastq.gz
│       ├── SRR1586885_pass_1.fastq.gz
│       ├── SRR1586886_pass_1.fastq.gz
│       ├── SRR1586887_pass_1.fastq.gz
│       ├── SRR1586888_pass_1.fastq.gz
│       ├── SRR1586889_pass_1.fastq.gz
│       ├── SRR1586890_pass_1.fastq.gz
│       ├── SRR1586891_pass_1.fastq.gz
│       ├── SRR1586892_pass_1.fastq.gz
│       ├── SRR1586893_pass_1.fastq.gz
│       ├── SRR1586894_pass_1.fastq.gz
│       ├── SRR1586895_pass_1.fastq.gz
│       ├── SRR1586896_pass_1.fastq.gz
│       ├── SRR1586897_pass_1.fastq.gz
│       ├── SRR1586898_pass_1.fastq.gz
│       ├── SRR1586899_pass_1.fastq.gz
│       ├── SRR1586900_pass_1.fastq.gz
│       ├── SRR1586901_pass_1.fastq.gz
│       ├── SRR1586902_pass_1.fastq.gz
│       ├── SRR1586903_pass_1.fastq.gz
│       ├── SRR1586904_pass_1.fastq.gz
│       ├── SRR1586905_pass_1.fastq.gz
│       ├── SRR1586906_pass_1.fastq.gz
│       ├── SRR1586907_pass_1.fastq.gz
│       ├── SRR1586908_pass_1.fastq.gz
│       ├── SRR1586909_pass_1.fastq.gz
│       ├── SRR1586910_pass_1.fastq.gz
│       ├── SRR1586911_pass_1.fastq.gz
│       ├── SRR1586912_pass_1.fastq.gz
│       ├── SRR1586913_pass_1.fastq.gz
│       ├── SRR1586914_pass_1.fastq.gz
│       ├── SRR1586915_pass_1.fastq.gz
│       └── SRR1586916_pass_1.fastq.gz
├── kallisto_info.txt
├── leaf_srr_numbers.txt
├── leaves_sra_info.txt
├── limmaDE.html
├── limmaDE.Rmd
├── logs
│   └── all.logs.tar.gz
├── merged_geno_sig.pdf
├── PRJNA262181_leaf_kallisto_info.txt
├── PRJNA262181_leaf_v3_kallisto_info.txt
├── README.md
├── Rplots.pdf
├── scripts
│   ├── kallisto_index.sh
│   ├── kallisto.sh
│   ├── kallisto_v3.sh
│   ├── single_kallisto.sh
│   └── sra-dump.sh
├── sig_genes_both.pdf
├── sig_genes_both.png
├── simple_meta.txt
├── sleuthAgg.Rmd
├── sra_info.txt
├── SraRunTable.txt
├── SRR1586620_info.txt
├── SRR1586766_info.txt
├── srr_numbers.txt
└── v3_b100_kallisto_info.txt

355 directories, 1224 files
