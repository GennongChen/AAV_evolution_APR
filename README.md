# AAV capsid evolution in human T cells identifies a family of CD7 targeting variants for efficient T and NK cell engineering
## Abstract
  This repository stores the code for the article “AAV capsid evolution in human T cells identifies a family of CD7 targeting variants for efficient T and NK cell engineering”, ensuring that other researchers can replicate all analyses and figures presented in the paper.
  We will provide data generated in our laboratory, public datasets, and analysis scripts. This includes four sections: AAV capsid sequence evolutionary analysis, CRISPR screening for AAV.APR31 receptor, RNA-seq for safety assessment, and CD7 distribution in human tissues.  
  **  To reproduce the code, please modify the tool and file paths in the scripts.**
## Data & Code
### Section 1: AAV capsid sequence evolutionary analysis (Fig.1 & Fig.S1)
  Raw sequencing and processed data are stored in GSE302084. And `qc_CapIV.sh` is used for preprocessing fastq files. `stat_CapIV.R` performs statistical analysis and visualization of AAV capsid sequences in the libraries using R. `20231219_CapIV_metadata.xlsx` contains the metadata corresponding to the raw fastq files.
### Section 2: CRISPR screening for AAV.APR31 receptor (Fig.4 & Fig.S3)
  Raw sequencing and processed data are stored in GSE302084.
### Section 3: RNA-seq for safety assessment (Fig.5 & Fig.S5)
  Raw sequencing data are stored in HRA012600. 
### Section 4: CD7 distribution in human tissues (Fig.6 & Fig.S6)
  HPA data was downloaded from HPA database (https://www.proteinatlas.org/humanproteome/tissue/data#hpa_tissues_rna). Cross-tissue multicellular single-cell atlas was downloaded from Zenodo (https://zenodo.org/records/15169362). ScRNA-seq data of bone marrow, thymus, and brain were downloaded from DISCO database (https://disco.bii.a-star.edu.sg/download). 

## Requirements
### tools
```
Fastp v0.234
Cutadapt v4.4.86
sickle v1.33
SeqKit 2.6.1
Mageck v0.5.9.5
STAR v2.5.2b
featureCounts v2.0.8
python v3.9.18
```
### python lib
```
anndata v0.9.2
pandas v2.1.1
scanpy v1.10.1
seaborn v0.13.2
```
### R lib
```
data.table v1.15.4
tidyverse v2.0.0
dplyr v1.1.4
ggplot2 v3.4.4
ggpubr v0.6.0
ggseqlogo v0.2
clusterProfiler v4.2.2
Pheatmap v1.0.12
factoextra v1.0.7
DESeq2 v1.42.0
Seurat4 v4.4.0
```
