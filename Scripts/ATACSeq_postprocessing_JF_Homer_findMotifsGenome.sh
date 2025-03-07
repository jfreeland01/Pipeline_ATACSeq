#!/bin/bash

### By: Jack Freeland (https://github.com/jfreeland01)
### Execute this to run HOMER findMotifsGenome.pl

parent_dir=""
input_dir="$parent_dir/DESeq2/homer" # wherever your input files are
nthread=""

### Loop over each .txt file in the input_dir
for file in "$input_dir"/*.txt; do
  # Extract the base name of the file (without the directory and extension)
  base_name=$(basename "$file" .txt)

  echo "Running Homer findMotifsGenome.pl on $file"

  # Define the output directory
  output_dir="$input_dir/$base_name"

  # Run the Homer command
  findMotifsGenome.pl "$file" hg38 "$output_dir" -size given -p "$nthread"

done