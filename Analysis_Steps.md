# Analysis Steps
## About
This document illustrates the individual steps of the **Computational Protocol for Assembly and Analysis of SARS-nCoV-2 Genomes**.
Please see the [README](https://github.com/banijolly/vslab-ncov2019-genome/blob/master/README.md) document for installation instructions. 


## Workflow

### Quality Assessment and Trimming
Tools used: FastQC, Trimmomatic.
```bash
fastqc <fastq file>
```
If the quality of the reads are not satisfactory, use Trimmomatic to remove bad quality reads and adapters.
For unpaired reads:
``` bash 
trimmomatic SE <input.fastq> <output.fastq>  ILLUMINACLIP:~anaconda3/envs/covid19-genome/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:100 
```
For paired reads:
``` bash 
trimmomatic PE <input_fwd.fastq> <input_rev.fastq> <output_fwd_paired.fastq> <output_fwd_unpaired.fastq> <output_rev_paired.fastq> <output_rev_unpaired.fastq> ILLUMINACLIP:~anaconda3/envs/covid19-genome/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:100 
```

### Analysis of species diversity present in the RNA sample
Tools used: Kraken, Krona
For unpaired reads:
```bash
kraken2 --db minikraken_8GB_20200312 --threads 16 --report <file.kreport> <trimmed.fastq> > <trimmed.kraken>
ktImportTaxonomy -s 3 -t 4 -o <visualization_output>.html <trimmed.kraken>
```
For paired reads:
``` bash kraken2 --db minikraken_8GB_20200312 --threads 16 --report <file.kreport> --paired <file1.fastq> <file2.fastq> > <trimmed.kraken>
ktImportTaxonomy -s 3 -t 4 -o <visualization_output>.html <trimmed.kraken>
```
<img src="https://github.com/banijolly/vslab-ncov2019-genome/blob/master/docs/Krona_output.png" align="center" height="100">




