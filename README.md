<p align="center">
  <img width="460" height="300" src="https://github.com/banijolly/vslab-ncov2019-genome/blob/master/docs/genepi_logo.png">
</p>

## About
The **Computational Protocol for Assembly and Analysis of SARS-nCoV-2 Genomes** has been compiled by VS-Lab at [CSIR-Insitute of Genomics an Integrative Biology](https://www.igib.res.in/) as an effort to aid analysis and interpretation of the sequencing data of SARS-CoV-2 using easy-to-use open source utilities using both reference-guided and de novo based strategies.
More information about out the lab and our work on COVID-19 can be found at the [lab website](http://vinodscaria.rnabiology.org/).

## Quickstart

### Installation
To use conda, download and install the [latest version of Anaconda](https://www.anaconda.com/distribution/).

Create and activate the covid19-genepi conda environment:
```bash
conda env create -f covid19-environment.yml
conda activate covid19-genepi
```
### Update Krona Taxonomy
Krona taxonomy databases will have to be manually updated before Krona can generate taxonomic report. The following code assumes Anaconda is installed in the home directory. The path can be updated according to your installation. 
```bash
bash ~/anaconda3/envs/covid19-genome/opt/krona/updateTaxonomy.sh 
bash ~/anaconda3/envs/covid19-genome/opt/krona/updateAccessions.sh
```

### Set up Minikraken database
The Minikraken database having complete bacterial, archaeal, and viral genomes in RefSeq is available for download at the [Kraken website](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads). 
```bash
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar -xvf minikraken_8GB_202003.tgz
export KRAKEN2_DB_PATH="<path/to/folder/containing/minikraken/database>"
```

### Install MEGAX
Install [MEGAX](https://www.megasoftware.net/) Command-Line Interface for analyze molecular evolution and generate phylogenetic trees:
```bash
wget https://www.megasoftware.net/do_force_download/megacc_10.1.1_amd64_beta.tar.gz
tar -zxvf megacc_10.1.1_amd64_beta.tar.gz
```
### Download Reference Genomes
The latest version of the human genome can be downloaded from [GENCODE]https://www.gencodegenes.org/human/ 
SARS-CoV2 genome can be downloaded from NCBI accession number [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)

### Download the test dataset
The RNA sequencing data of 14 patients infected with SARS-CoV-2 sequenced by University of Washington can be download from [SRA repository](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP251618).

### Download Reference Dataset for Phylogenetic Analysis
The [Global Initiative on Sharing All Influenza Data (GISAID)](https://www.gisaid.org/) gives public access to the most complete repository of sequencing data for SARS-CoV2. The sequences for phylogenetic analysis can be downloaded from the [EpiCoV portal](https://www.epicov.org/epi3/) after creating an account and signing in.
