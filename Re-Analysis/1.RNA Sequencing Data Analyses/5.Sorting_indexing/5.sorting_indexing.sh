#!/bin/bash

#SBATCH --job-name=samtools
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=20:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ra500@student.le.ac.uk

module load samtools/1.17-wenuvv5

# Input directory containing the .sam files
input_dir="/scratch/alice/r/ra500/alignement/"

# Output directory for the sorted BAM files
output_dir="/scratch/alice/r/ra500/sorted_indexed"

# Navigate to the input directory
cd "$input_dir" || exit

# Loop over each .sam file
for sam_file in *.sam; do
    	# Extract the base name of the file (without extension)
    	base_name=$(basename "$sam_file" .sam)
   	 
    	# Sort the .sam file using samtools and save to the output directory
    	samtools sort "$sam_file" -o "${output_dir}/${base_name}_sorted.bam"
   	 
    	# Optional: Index the sorted BAM file
    	samtools index "${output_dir}/${base_name}_sorted.bam"
done
