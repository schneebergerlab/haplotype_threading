#!/bin/bash

# Check if a file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <file.txt>"
    exit 1
fi

# Input file
input_file=$1

# Extract the base name without extension
base_name=$(basename "$input_file" .txt)
# Output file
output_file="${base_name}.fasta"

# Check if the file exists
if [ ! -f "$input_file" ]; then
    echo "File not found!"
    exit 1
fi

# Initialize a counter for unique IDs
counter=1

#Bash
awk -F'\t' '{
    printf ">%d_%s\n%s\n", NR, $1, $2
}' "$input_file" > "$output_file"

echo "FASTA file created: $output_file"

