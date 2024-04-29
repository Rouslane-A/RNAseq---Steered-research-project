#!/bin/bash

# Set the directory containing the .txt files
folder="/scratch/alice/r/ra500/fastq_screen_out/"

# Create an associative array to keep track of processed prefixes
declare -A processed_prefixes

# Loop through each .txt file in the folder
for file in "$folder"/*.txt; do
    # Extract filename without extension
    filename=$(basename "$file")
    filename_no_ext="${filename%.*}"

    # Extract the prefix before the first '_'
    prefix="${filename_no_ext%%_*}"

    # Check if the prefix has already been processed
    if [[ ! -v processed_prefixes[$prefix] ]]; then
        # Process the file
        awk -v fname="$filename" '
            BEGIN {
                FS="\t"  # Set the field separator to tab
            }

            NR==3 {
                reads = $2
                one_hit_one_genome_h = $4
                multiple_hits_one_genome_h = $5
                one_hit_multiple_genomes_h = $6
                multiple_hits_multiple_genomes_h = $7
            }
            
            NR==4 {
                one_hit_one_genome_m = $4
                multiple_hits_one_genome_m = $5
                one_hit_multiple_genomes_m = $6
                multiple_hits_multiple_genomes_m = $7
            }

            END {
                if ((one_hit_one_genome_h + multiple_hits_one_genome_h) >= 5 * (one_hit_one_genome_m + multiple_hits_one_genome_m)) {
                    split(fname, parts, /_/)
                    print parts[1]
                }
            }
        ' "$file"

        # Mark the prefix as processed
        processed_prefixes[$prefix]=1
    fi
done > output.txt


