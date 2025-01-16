#!/bin/bash

### By: Jack Freeland (https://github.com/jfreeland01)
### Execute this to get consensus peak files files

### Directories 
parent_dir=""
peak_dir="$parent_dir"/macs3
bam_dir="$parent_dir"/bam

##### Generate concensus peak file #####

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] -------------- Merging Peaks --------------"

### Find all narrowPeak files and concatenate them into one file
find "$peak_dir" -name "*.narrowPeak" -exec cat {} + > "$peak_dir/all_concatenate.narrowPeak.bed"

#### Sort the combined file
sort -k1,1 -k2,2n "$peak_dir/all_concatenate.narrowPeak.bed" > "$peak_dir/all_concatenate_sorted.narrowPeak.bed"

#### Merge the peaks
bedtools merge -i "$peak_dir/all_concatenate_sorted.narrowPeak.bed" > "$peak_dir/all_concatenate_sorted_merged.narrowPeak.bed"

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] -------------- Finished Merging Peaks --------------"


##### Generate counts file #####

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] -------------- Counting Reads in Peaks --------------"

# Create an array of all BAM files
bam_files=($bam_dir/*V8.bam)


# Create header with BAM filenames for the multicov output, ensuring tab separation
header="chr\tstart\tend"
for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" | sed 's/_V8.bam//g')
    header="$header\t$sample_name"
done

# Output the header to the file
echo -e "$header" > "$peak_dir/multicov_merged_counts.txt"

# Run bedtools multicov to count reads in merged peaks across all BAM files
bedtools multicov -bams "${bam_files[@]}" -bed "$peak_dir/all_concatenate_sorted_merged.narrowPeak.bed" >> "$peak_dir/multicov_merged_counts.txt"

timestamp=$(date '+%Y-%m-%d %H:%M:%S')
echo "[$timestamp] -------------- Finished Counting Reads --------------"