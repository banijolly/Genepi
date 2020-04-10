# Analysis Steps
## About
This document illustrates the individual steps of the **Computational Protocol for Assembly and Analysis of SARS-nCoV-2 Genomes**.
Please see the [README](https://github.com/banijolly/vslab-ncov2019-genome/blob/master/README.md) document for installation instructions. 


## Workflow

### Quality Assessment and Trimming
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
For unpaired reads:
```bash
kraken2 --db minikraken_8GB_20200312 --threads 16 --report <file.kreport> <trimmed.fastq> > <trimmed.kraken>
ktImportTaxonomy -s 3 -t 4 -o <visualization_output>.html <trimmed.kraken>
```
For paired reads:
``` bash kraken2 --db minikraken_8GB_20200312 --threads 16 --report <file.kreport> --paired <file1.fastq> <file2.fastq> > <trimmed.kraken>
ktImportTaxonomy -s 3 -t 4 -o <visualization_output>.html <trimmed.kraken>
```
<img src="https://github.com/banijolly/vslab-ncov2019-genome/blob/master/docs/Krona_output.png" height="250">

### Reference based assembly against the human reference genome
For unpaired reads,
```bash 
hisat2 -x <Humangenomeprefix> -U <unpaired_input.fastq> -S <human.sam> -p 16 --dta-cufflinks --summary-file <humanUnpaired.log>
```

For paired reads,
```bash
hisat2 -x <Humangenomeprefix> -1 <file1.fastq> -2 <file2.fastq> -S <human.sam> -p 16 --dta-cufflinks --summary-file <file_human_paired.log>
```

###Pre-processing and extracting unaligned reads
```bash
samtools sort -@ 8 <human.sam> -o <human.bam> -O BAM
samtools view -u -f 4 -F264 <human.bam> > temp1.bam
samtools view -u -f 8 -F 260 <human.bam> > temp2.bam
samtools view -u -f 12 -F 256 <human.bam> > temp3.bam
samtools merge -u - temp[123].bam | samtools sort -n -o <unmapped.bam>
```
####Converting BAM to FASTQ
For unpaired reads,
```bash
bamToFastq -i <unmapped.bam> -fq <unmapped.fq>
```

For paired reads,
```bash
bamToFastq -i <unmapped.bam> -fq <unmapped_1.fq> -fq2 <unmapped_2.fq>
```

###Alignment to the viral reference genome

For unpaired reads,
```bash
hisat2 -x <NC_045512.2_viral_reference_genome> -U <unmapped.fq> -S <covid.sam> -p 16 --dta-cufflinks --summary-file <unmapped.log>
```

For paired reads,
```bash hisat2 -x <NC_045512.2_viral_reference_genome> -1 <unmapped_1.fq> -2 <unmapped_2.fq> -S <covid.sam> -p 16 --dta-cufflinks --summary-file <unmapped.log>
```

####Evaluation of Alignment statistics to the viral reference genome

```bash 
samtools flagstat <covid.sam>
samtools sort -@ 8 <covid.sam> -o <covid.bam> -O BAM
```



