#!/bin/bash

### By: Jack Freeland (jackfreeland01@gmail.com, https://www.linkedin.com/in/jack-freeland-384526142)
### Execute this to run computeMatrix & plotHeatmap
### File names (file1/2/3) and what you are search for (i.e. *.wig) will need to be modified

parent_dir=""
wig_dir="$parent_dir/Wig"
ref_dir="$parent_dir/references"
bigwig_dir="$parent_dir/bigWig"

file1=""
file2=""
file3=""

### Generate mean wig file

mkdir -p "$wig_dir/mean"

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] Starting wiggletools mean for $file1 & $file2 & $file3"

wiggletools mean "$wig_dir/${file1}.wig" "$wig_dir/${file2}.wig" "$wig_dir/${file3}.wig" > "$wig_dir/${file1}_${file2}_${file3}_mean.wig"

timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
echo "[$timestamp] Finished wiggletools mean"


### Convert wig to bigwig file

# generated reference file from index fasta file by the following 
# cut -f1,2 [fasta-file].fai > grch38_noalt_decoy_as.chrom.sizes

mkdir -p "$bigwig_dir/mean"

for wig_file in "$wig_dir"/*.wig; do

	sample_ID=$(basename "$wig_file" .wig)

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] Converting from wig to bigwig: $sample_ID"

    wigToBigWig "$wig_file" "$ref_dir/GRCH38_noalt_decoy_as.chrom.sizes" "$bigwig_dir/${sample_ID}.bw"

done