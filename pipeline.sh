#!/bin/bash 
start=`date +%s`
usage()
{
  echo "Usage: ./run_PE.sh -d <Directory of fastq files>  -r <Covid_reference_genome.fasta> -z <suffix of samples (_R1.fastq.gz/_R2.fastq.gz or _R1.fastq/_R2.fastq)  -o <Output directory suffix> -h <help>"
  exit 2
}

while getopts d:r:z:o:h: option 
do 
 case "${option}" 
 in 
 d) DIRECTORY=${OPTARG};;
 r) COVID_REFERENCE=${OPTARG};;
 z) FASTQ_ZIP=${OPTARG};;
 o) OUT=${OPTARG};;
 h|?) usage ;; esac
done

mkdir Fastqc_output_$OUT
mkdir Trimmed_output_$OUT
mkdir Hisat2_covid_output_$OUT
mkdir Variant_calling_output_$OUT
mkdir Picard_metrics
mkdir PANGO_reports 
touch SAMPLE_IDs

printf "Sample ID\tLineage\tTotal Reads R1\tTotal Reads R2\tTrimmed Reads R1\tTrimmed Reads R2\tHISAT2_alignment_percentage\tX coverage\tGenome coverage\tVCF Variant Count\tNs in FASTA\n" > SampleSummary.txt

FileList="$(ls $DIRECTORY/*$FASTQ_ZIP | awk '{ print $1}' | awk -F'/' '{ print $2}' | grep "_R1" | awk -F'_R1' '{ print $1}')"


for i in $FileList
do
	id=$(echo $i | sed 's/_S[0-9].*//g')
	echo $id >> SAMPLE_IDs
	echo ------------------------------------------------------------------
	echo "Starting analysis for "$i
	echo ------------------------------------------------------------------

	echo -e "Scanning Location:"$DIRECTORY"\t"$i
	echo ------------------------------------------------------------------
	echo "Evaluating FASTQC report for Sample"
	echo ------------------------------------------------------------------

	FASTQC_COMMAND_R1="fastqc $DIRECTORY/"$i"_R1_"$FASTQ_ZIP" -o Fastqc_output_"$OUT
	echo $FASTQC_COMMAND_R1
	eval "$FASTQC_COMMAND_R1"
	FASTQC_COMMAND_R2="fastqc $DIRECTORY/"$i"_R2_"$FASTQ_ZIP" -o Fastqc_output_"$OUT
	echo $FASTQC_COMMAND_R2
	eval "$FASTQC_COMMAND_R2"
	echo ------------------------------------------------------------------
	echo "FASTQC Evaluation Completed"
	echo ------------------------------------------------------------------



	echo ------------------------------------------------------------------
	echo "Using Trimmomatic to trim low quality bases"
	echo ------------------------------------------------------------------

	TRIM_COMMAND_PE="trimmomatic PE $DIRECTORY"/""$i"_R1_"$FASTQ_ZIP" $DIRECTORY"/""$i"_R2_"$FASTQ_ZIP" Trimmed_output_"$OUT"/"$i"_fwd_paired.fastq Trimmed_output_"$OUT"/"$i"_fwd_unpaired.fastq Trimmed_output_"$OUT"/"$i"_rev_paired.fastq Trimmed_output_"$OUT"/"$i"_rev_unpaired.fastq ILLUMINACLIP:"$HOME"/anaconda3/envs/covid19-genepi/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:35"	
echo $TRIM_COMMAND_PE
	eval "$TRIM_COMMAND_PE"

	echo ------------------------------------------------------------------
	echo "Trimmomatic run completed"
	echo ------------------------------------------------------------------


	echo ------------------------------------------------------------------
	echo "Building HISAT2 indexes for reference genome"
	echo ------------------------------------------------------------------

	HISAT_BUILD_COMMAND="hisat2-build "$COVID_REFERENCE" "$COVID_REFERENCE
	echo $HISAT_BUILD_COMMAND
	eval "$HISAT_BUILD_COMMAND"

	echo ------------------------------------------------------------------
	echo "Done Building Indexes for "$COVID_REFERENCE
	echo ------------------------------------------------------------------
        
	echo ------------------------------------------------------------------
	echo "Using HISAT2 to perform reference based mapping"
	echo ------------------------------------------------------------------

	HISAT_COVID_COMMAND="hisat2 -x "$COVID_REFERENCE" -1 Trimmed_output_"$OUT"/"$i"_fwd_paired.fastq -2 Trimmed_output_"$OUT"/"$i"_rev_paired.fastq -S Hisat2_covid_output_"$OUT"/"$i"_aligned.sam -p 16 --summary-file Hisat2_covid_output_"$OUT"/"$i"_hisat.log"        
	echo $HISAT_COVID_COMMAND
	eval "$HISAT_COVID_COMMAND"
	echo ------------------------------------------------------------------
	echo "Reference mapping using HISAT2 completed"
	echo ------------------------------------------------------------------



	echo ------------------------------------------------------------------
	echo "Converting SAM to BAM and performing coordinate sorting"
	echo ------------------------------------------------------------------

	SAM2BAM_COMMAND="samtools sort Hisat2_covid_output_"$OUT"/"$i"_aligned.sam -o Hisat2_covid_output_"$OUT"/"$i"_aligned.bam"
	echo "$SAM2BAM_COMMAND"
	eval "$SAM2BAM_COMMAND"
	echo ------------------------------------------------------------------
	echo "Post-processing completed"
	echo ------------------------------------------------------------------

	COVID_FLAGSTAT_COMMAND="samtools flagstat Hisat2_covid_output_"$OUT"/"$i"_aligned.bam > Hisat2_covid_output_"$OUT"/"$i"_flagstat.txt"
	echo "$COVID_FLAGSTAT_COMMAND"
	eval "$COVID_FLAGSTAT_COMMAND"

	echo ------------------------------------------------------------------
	echo "Generating pileup file for variant calling"
	echo ------------------------------------------------------------------

	COVID_MPILEUP="samtools mpileup -f "$COVID_REFERENCE" Hisat2_covid_output_"$OUT"/"$i"_aligned.bam > Variant_calling_output_"$OUT"/"$i".pileup"
	echo "$COVID_MPILEUP"
	eval "$COVID_MPILEUP"
	echo ------------------------------------------------------------------
	echo "Pileup file generated"
	echo ------------------------------------------------------------------

	echo ------------------------------------------------------------------
	echo "Calling Variants using VarScan"
	echo ------------------------------------------------------------------

	VARSCAN_COMMAND="varscan mpileup2cns Variant_calling_output_"$OUT"/"$i".pileup --output-vcf 1 --variants > Variant_calling_output_"$OUT"/"$i".vcf"
	echo "$VARSCAN_COMMAND"
	eval "$VARSCAN_COMMAND"
	echo ------------------------------------------------------------------
	echo "Variant calling using VarScan Completed"
	echo ------------------------------------------------------------------

	echo ------------------------------------------------------------------
	echo "Generating FASTA Sequence for the sample"
	echo ------------------------------------------------------------------

	CONSENSUS_MPILEUP_COMMAND="samtools mpileup -uf "$COVID_REFERENCE" Hisat2_covid_output_"$OUT"/"$i"_aligned.bam | bcftools call -c | vcfutils.pl vcf2fq > Variant_calling_output_"$OUT"/"$i"_consensus.fq"
	eval "$CONSENSUS_MPILEUP_COMMAND"

	FASTQ_FASTA_COMMAND="seqtk seq -aQ64 -q20 -n N Variant_calling_output_"$OUT"/"$i"_consensus.fq > Variant_calling_output_"$OUT"/"$i"_consensus.fasta"
	echo "$FASTQ_FASTA_COMMAND"
	eval "$FASTQ_FASTA_COMMAND"
	echo ------------------------------------------------------------------
	echo "Generated FASTA Sequence for the sample. The output is stored in : Variant_calling_output_"$OUT"/"$i"_consensus.fasta"
	echo ------------------------------------------------------------------

	echo ------------------------------------------------------------------
	echo "Running picard to collect coverage information"
	echo ------------------------------------------------------------------
	PICARD_COMMAND="picard CollectMultipleMetrics -I Hisat2_covid_output_"$OUT"/"$i"_aligned.bam -O Picard_metrics/"$i" -R "$COVID_REFERENCE
	eval "$PICARD_COMMAND"
	echo ------------------------------------------------------------------
	echo "Picard run completed"
	echo ------------------------------------------------------------------

	echo ------------------------------------------------------------------
	echo "Generating Report for Sample:"$i
	echo ------------------------------------------------------------------

	r1_cmd="zcat "$DIRECTORY"/"$i"_R1_"$FASTQ_ZIP" | wc -l | awk -F'\t' '(\$1=\$1/4)'"
	r1=$(eval "$r1_cmd")

	r2_cmd="zcat "$DIRECTORY"/"$i"_R2_"$FASTQ_ZIP" | wc -l | awk -F'\t' '(\$1=\$1/4)'"
	r2=$(eval "$r2_cmd")

	tr1_cmd="cat Trimmed_output_"$OUT"/"$i"_fwd_paired.fastq | wc -l |  awk -F'\t' '(\$1=\$1/4)'"
	tr1=$(eval "$tr1_cmd")

	tr2_cmd="cat Trimmed_output_"$OUT"/"$i"_rev_paired.fastq | wc -l |  awk -F'\t' '(\$1=\$1/4)'"
	tr2=$(eval "$tr2_cmd")

	hsa_cmd="tail -1 Hisat2_covid_output_"$OUT"/"$i"_hisat.log | sed 's/\% overall alignment rate//g'"  
	hsa=$(eval "$hsa_cmd")

	bpcov_cmd="cat Picard_metrics/"$i".alignment_summary_metrics  | grep -w "PAIR" | cut -f8"
	bpcov=$(eval "$bpcov_cmd")

	xcov=$(echo $bpcov| awk '($1=$1/29903)')

	grzero_cmd="samtools depth -a Hisat2_covid_output_"$OUT"/"$i"_aligned.bam | awk -F'\t' '(\$3>0)' | wc -l"
	grzero=$(eval "$grzero_cmd")

	genomecov=$(echo $grzero | awk '($1=100*$1/29903)')

	varcount_cmd="cat Variant_calling_output_"$OUT"/"$i".vcf  | grep -v '#' | wc -l "
	varcount=$(eval "$varcount_cmd")

	fastaN_cmd="cat Variant_calling_output_"$OUT"/"$i"_consensus.fasta  | grep -v '>' | grep -io N| wc -l "
	fastaN=$(eval "$fastaN_cmd")

	echo ------------------------------------------------------------------
	echo "Assigning lineage using Pangolin"
	echo ------------------------------------------------------------------

	pangolin "Variant_calling_output_"$OUT"/"$i"_consensus.fasta" --outdir PANGO_reports --outfile "$i"_lineage-report.csv
	lineage_cmd="cat PANGO_reports/"$i"_lineage-report.csv | cut -d ',' -f2 | tail -1"
	lineage=$(eval "$lineage_cmd")
	

	printf $id"\t$lineage\t"$r1"\t"$r2"\t"$tr1"\t"$tr2"\t"$hsa"\t"$xcov"\t"$genomecov"\t"$varcount"\t"$fastaN"\n" >> SampleSummary.txt
	echo ------------------------------------------------------------------
	echo "Analysis completed for Sample:"$i
	echo ------------------------------------------------------------------

done

echo ------------------------------------------------------------------
echo "Combining all FASTA Files and storing output in Combined_Fasta_sequences.fasta"
echo ------------------------------------------------------------------


mkdir FASTAS

for i in `cat SAMPLE_IDs`
do
cp Variant_calling_output_*/"$i"_*fasta FASTAS/"$i"_consensus.fasta
done

cd FASTAS
python ../change_header.py .
rm -rf *fasta
cat *fa > Combined_Fasta_sequences.fasta

rm -rf *fa



sed 's/_consensus//g' Combined_Fasta_sequences.fasta | fold -100 > temp
mv temp ../Combined_Fasta_sequences.fasta
cd ../
rm -rf FASTAS



end=`date +%s`

runtime=$((end-start))

echo ------------------------------------------------------------------
echo All done
echo Time Taken: $runtime seconds
echo ------------------------------------------------------------------


