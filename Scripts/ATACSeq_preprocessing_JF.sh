#!/bin/bash

### Pipeline for preprocessing ATAC-Seq data
### By: Jack Freeland (https://github.com/jfreeland01)
### This pipeline will take take a directory ($parent_dir) which contains a folder of raw fastq files ($raw_dir) from a paired-end ATACSeq experiment and generate BAM/wig/bigwig files and call peaks.
### *NOTE* - Depending on your fastq file naming format, the method used to generate $sample_ID may need to be modified.

### Provide a path for the following directories:

### $parent_dir - directory for all output files
### $raw_dir - sub-directory containing raw fastq files
### $ref_dir - sub-directory containing genome reference/index for bowtie2 (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
### $nthread - number of threads to be run. Note - When filtering for PCR dups, it asks for number of threads (default nthread) and an amount of memory to be alocated to each (default 3G Xmx3G). So ensure you have 3Gb*nthread of memory available. If not, lower memory allocation (e.g., 1-2Gb) or replace nthread with something lower value for function (lowering memory and keeping cores higher seems to perform better)


parent_dir=""
raw_dir="$parent_dir/raw_fastq"
ref_dir="$parent_dir/references"
nthread=""


##### FASTQC ON RAW READS #####
fastqc_raw_dir="$parent_dir/fastqc/raw"

for subdir in "$raw_dir"/*; do

	sample_ID=$(basename "$subdir" | cut -d '_' -f 1)
	mkdir -p "$fastqc_raw_dir/$sample_ID"

	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Running fastqc on raw $sample_ID"

    fastqc -o "$fastqc_raw_dir/$sample_ID" "$subdir"/*fastq* \
	-t "$nthread"

done


##### CUTADAPT ON RAW READS (TRIM ADAPTER) #####
cutadapt_dir="$parent_dir/cutadapt_fastq"

# adapter sequence varies depending on kit used for library prep. Some wrapprs of cutadapt (e.g., TrimGalore) can automatically detect which is which (they also list some of the common ones)
# https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

for subdir in "$raw_dir"/*; do

	sample_ID=$(basename "$subdir" | cut -d '_' -f 1)
	mkdir -p "$cutadapt_dir/$sample_ID"
	
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Running cutadapt on raw $sample_ID"

    python3 -m cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA \
	-j "$nthread" -q 20 \
	-O 6 \
	-m 35 \
	-o "$cutadapt_dir/$sample_ID/${sample_ID}_R1.fastq.gz" \
	-p "$cutadapt_dir/$sample_ID/${sample_ID}_R2.fastq.gz" \
	"$subdir"/*_R1_* "$subdir"/*_R2_*

done


##### FASTQC ON TRIMMED READS #####
fastqc_cutadapt_dir="$parent_dir/fastqc/cutadapt"

for subdir in "$cutadapt_dir"/*; do

	sample_ID=$(basename "$subdir")
	mkdir -p "$fastqc_cutadapt_dir/$sample_ID"

	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Running fastqc on trimmed $sample_ID"
        
	fastqc -o "$fastqc_cutadapt_dir/$sample_ID" "$subdir"/*fastq* \
	-t "$nthread"

done


##### BOWTIE2 ON TRIMMED READS (GENERATE .SAM FILES) #####
sam_dir="$parent_dir/sam"

# ref: https://benlangmead.github.io/aws-indexes/bowtie

mkdir -p $sam_dir

for subdir in "$cutadapt_dir"/*; do
  
	sample_ID=$(basename "$subdir")
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Running Bowtie2 on $sample_ID"

  bowtie2 --non-deterministic --mm --phred33 --very-sensitive \
  -p "$nthread" \
  -x "$ref_dir/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as" \
  -1 "$subdir"/*R1.fastq* \
  -2 "$subdir"/*R2.fastq* \
  -S "$sam_dir/${sample_ID}.sam"
done


##### CONVERT .SAM TO .BAM FILES #####
bam_dir="$parent_dir/bam"

mkdir -p "$bam_dir/raw"

for sam_file in "$sam_dir"/*.sam; do

	sample_ID=$(basename "$sam_file" .sam)
    
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Converting ${sample_ID}.sam to ${sample_ID}.bam"
	
	samtools view -@ "$nthread" -bS "$sam_file" > "$bam_dir/raw/${sample_ID}.bam"

	rm "$sam_file"

done


##### FILTER .BAM FILES #####
mkdir -p "$bam_dir/processed"
mkdir -p "$bam_dir/dup_files"

for bam_file in "$bam_dir"/raw/*.bam; do

	sample_ID=$(basename "$bam_file" .bam)

	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] -------------- STARTING $sample_ID --------------"
	
	### Sort reads
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Sorting reads: $sample_ID"
	
    if [ ! -e "$bam_file" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to sort: $sample_ID"
        continue  # Skip this iteration and proceed to the next $bam_file. Intended to stop the loop from rm previous bam file without generating the next 'version'
    fi

	samtools sort -@ "$nthread" -o "$bam_dir/processed/${sample_ID}_V1.bam" "$bam_file"
	bam_V1="$bam_dir/processed/${sample_ID}_V1.bam"
	
	### Index file
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Indexing: $sample_ID"
	
   if [ ! -e "$bam_V1" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to index: $bam_V1"
        continue
    fi

	samtools index "$bam_V1"

	### Adjust for transposase binding offset
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Adjust for transposase binding offset: $sample_ID"

	alignmentSieve --bam "$bam_V1" \
	-o "$bam_dir/processed/${sample_ID}_V2.bam" \
    --ATACshift --numberOfProcessors "$nthread"

	rm "$bam_V1"
	rm "${bam_V1}.bai"
	bam_V2="$bam_dir/processed/${sample_ID}_V2.bam"

	### Filter unaligned reads
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Filtering unaligned reads: $sample_ID"
	
    if [ ! -e "$bam_V2" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to filter unaligned reads: $bam_V2"
        continue
    fi

	samtools view -@ "$nthread" -b -F 4 -o \
	"$bam_dir/processed/${sample_ID}_V3.bam" "$bam_V2"

	rm "$bam_V2"
	bam_V3="$bam_dir/processed/${sample_ID}_V3.bam"
	
	### Filter ChrM and 'decoy' reads
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Filtering ChrM reads: $sample_ID"
	
    if [ ! -e "$bam_V3" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to filter ChrM reads: $bam_V3"
        continue
    fi

	samtools view -@ "$nthread" -h "$bam_V3" | \
	egrep -v "chrM|Un|random|decoy" | \
	samtools view -b -o "$bam_dir/processed/${sample_ID}_V4.bam"

	rm "$bam_V3"
	bam_V4="$bam_dir/processed/${sample_ID}_V4.bam"

	### Sort reads V2
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Sorting reads V2: $sample_ID"
	
    if [ ! -e "$bam_V4" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to sort: $bam_V4"
        continue 
    fi

	samtools sort -@ "$nthread" -o "$bam_dir/processed/${sample_ID}_V5.bam" "$bam_V4"
	rm "$bam_V4"
	bam_V5="$bam_dir/processed/${sample_ID}_V5.bam"


	### Filter PCR Duplicates
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Filtering PCR Duplicates: $sample_ID"

    if [ ! -e "$bam_V5" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to filter PCR dups: $sample_ID"
        continue
    fi

	java -Xmx3G -XX:ParallelGCThreads="$nthread" \
	-jar "$ref_dir/picard.jar" MarkDuplicates \
	-M "$bam_dir/dup_files/${sample_ID}_m.txt" \
	-I "$bam_V5" \
	-O "$bam_dir/processed/${sample_ID}_V6.bam" \
	-REMOVE_DUPLICATES true

	rm "$bam_V5"
	bam_V6="$bam_dir/processed/${sample_ID}_V6.bam"

	### Filter for quality reads
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Filtering Quality Reads >30: $sample_ID"

    if [ ! -e "$bam_V6" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to filter quality reads: $bam_V6"
        continue
    fi

	samtools view -b -q 30 -@ "$nthread" "$bam_V6" > "$bam_dir/processed/${sample_ID}_V7.bam"

	rm "$bam_V6"
	bam_V7="$bam_dir/processed/${sample_ID}_V7.bam"

	### Filter of ENCODE blacklists
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Filtering ENCODE blacklist regions: $sample_ID"

    if [ ! -e "$bam_V7" ]; then
        echo "[$timestamp] ERROR: Input BAM file not found to filter blacklist regions: $bam_V7"
        continue
    fi

	bedtools subtract -a "$bam_V7" -b "$ref_dir/hg38-blacklist.v2_sorted.bed" \
	-A > "$bam_dir/processed/${sample_ID}_V8.bam"

	rm "$bam_V7"

	echo "[$timestamp] -------------- FINISHED $sample_ID --------------"

done


##### INDEX .BAM FILES #####
for bam_file in "$bam_dir"/processed/*V8.bam; do

	sample_ID=$(basename "$bam_file" | cut -d '_' -f 1)
    
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Indexing $sample_ID"
	
	samtools index $bam_file

done


##### Generate bigWig files #####

# reference from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
bigWig_dir="$parent_dir/bigWig"
mkdir -p "$bigWig_dir"

for bam_file in "$bam_dir"/processed/*V8.bam; do

	sample_ID=$(basename "$bam_file" | cut -d '_' -f 1)
    
	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Generating bigWig file for $sample_ID"

	bamCoverage -b "$bam_file" --normalizeUsing RPGC \
	--effectiveGenomeSize 2913022398 \
	-p "$nthread" \
	-o "$bigWig_dir/${sample_ID}.bw" \
	--binSize 10

	timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Finished generating bigWig file for $sample_ID"

done


##### Generate Wig files from bigWig #####
Wig_dir="$parent_dir/Wig"
mkdir -p "$Wig_dir"

for bigwig_file in "$bigWig_dir"/*.bw; do

	sample_ID=$(basename "$bigwig_file" .bw)

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] Converting bigWig to Wig for $sample_ID"

    bigWigToWig "$bigWig_dir/${sample_ID}.bw" "$Wig_dir/${sample_ID}.wig"

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] Finished converting bigWig to Wig for $sample_ID"

done


##### Generate Peak files using Macs3 (parallel) #####
peak_dir="$parent_dir/macs3"

export parent_dir
export bam_dir
export peak_dir
export macs3

process_sample() {
    bam_file="$1"

    # Check if the filename is not empty
    if [[ -z "$bam_file" || ! -f "$bam_file" ]]; then
        echo "Finished with all bams in given directory"
        return 1  # Exit the function with an error status
    fi

    sample_ID=$(basename "$bam_file" | cut -d '_' -f 1)
    mkdir -p "$peak_dir/$sample_ID"

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] Running Macs3 on $sample_ID"

    # Run Macs3 callpeak
    macs3 callpeak -f BAMPE -g hs -q 0.01 -t "$bam_file" -n "$sample_ID" --outdir "$peak_dir/$sample_ID"
}

export -f process_sample

# Find all BAM files and pipe them into parallel
find "$bam_dir"/processed -name '*V8.bam' | parallel -j "$nthread" process_sample

# Ouput format Chromosome | Start Position | End Position | Name of the Peak | Score | Strand | Signal Value | p-value (-log10) | q-value (FDR) | Peak Point