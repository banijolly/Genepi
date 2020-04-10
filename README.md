# vslab-ncov2019-genome
## About
The **Computational Protocol for Assembly and Analysis of SARS-nCoV-2 Genomes** has been compiled by VS-Lab at [CSIR-Insitute of Genomics an Integrative Biology](https://www.igib.res.in/) as an effort to aid analysis and interpretation of the sequencing data of SARS-CoV-2 using easy-to-use open source utilities using both reference-guided and de novo based strategies.
More information about out the lab and our work on COVID-19 can be found at the [lab website](http://vinodscaria.rnabiology.org/).

## Quickstart

### Installation
To use conda, download and install the [latest version of Anaconda](https://www.anaconda.com/distribution/).

Create and activate the covid19-genome conda environment:
```bash
conda env create -f covid19-environment.yml
conda activate covid19-genome
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



