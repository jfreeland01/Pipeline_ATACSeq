#!/bin/bash

### By: Jack Freeland (https://github.com/jfreeland01)
### Execute this to run computeMatrix & plotHeatmap (block out chunks as all will not be needed at same time)

### Directories 
parent_dir=""
ref_dir="$parent_dir"/references
bigwig_dir="$parent_dir"/bigWig
matrix_dir="$parent_dir"/computeMatrix
nthread=""

mkdir -p "$matrix_dir"


### Make one matrix for all samples
timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] Running ComputeMatrix"

computeMatrix reference-point \
    --referencePoint TSS \
    --regionsFileName "$ref_dir"/TSS_1_V2.bed \
    --scoreFileName "$bigwig_dir"/* \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --binSize 10 \
    --outFileName "$matrix_dir/ComputeMatrix_RefPnt.gz" \
    --numberOfProcessors "$nthread" \
    --missingDataAsZero \
    --skipZeros

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] Finished Running ComputeMatrix"


### Make one matrix per sample (less likely relavent)
for bw_file in "$merged_bigwig_dir"/*; do

    sample_ID=$(basename "$bw_file" | cut -d '_' -f 1)

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Computing Matrix for $sample_ID"

    computeMatrix reference-point \
        --referencePoint TSS \
        --regionsFileName "$ref_dir"/TSS_1_V2.bed \
        --scoreFileName "$bw_file" \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --binSize 10 \
        --outFileName "$matrix_dir/${sample_ID}_ComputeMatrix_RefPnt.gz" \
        --numberOfProcessors "$nthread" \
        --missingDataAsZero \
        --skipZeros

done


### Plot heatmap
plotHeatmap -m "$matrix_dir"/ComputeMatrix_RefPnt_ATAC_.gz \
    -out "$matrix_dir"/OverallTSS_ATAC_heatmap.png



# To get this TSS file, I downloaded it from USCC genome browser

# Go to the UCSC Table Browser:

# This browser allows you to query and download the specific genomic data you need.
# Configure Your Query:

# clade: Mammal
# genome: Human (or choose your organism of interest)
# assembly: GRCh38/latest (or the genome assembly relevant to your study)
# group: Genes and Gene Predictions
# track: GENCODE v38 (or choose the version appropriate to your assembly)
# table: GENCODE v38 (this table contains gene annotations)
# region: genome (to get data for the entire genome)
# output format: selected fields from primary and related tables
# output file: Enter a name, e.g., TSS.bed
# Click get output.
# Select Fields for Output:

# After clicking get output, you'll be prompted to select the fields to output:
# Select chrom, txStart, txEnd, strand, name.
# Click get output again.

# Have to edit file to put in correct order

# awk 'BEGIN{FS="\t"; OFS="\t"} 
#      NR > 1 {if ($3 == "+") print $2, $4-1, $4, $1; 
#              else if ($3 == "-") print $2, $5-1, $5, $1}' TSS_1.bed > TSS_1_V2.bed