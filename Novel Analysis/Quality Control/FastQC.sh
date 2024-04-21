#!/bin/bash

# ALICE settings
#SBATCH --job-name=STAR_index
#SBATCH --nodes=2
#SBATCH --tasks-per-node=4
#SBATCH --mem=60gb
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bg171@student.le.ac.uk

# Load the required module
module load fastqc/0.12.1-hkgpcde

# Define input and output directories
input_dir="/scratch/alice/b/bg171/SRP_Data" # this is the sequenced cells we downloaded from the database
output_dir="/scratch/alice/b/bg171/rawreads_QA_stats" # this is the output of the QC data

# Change to the input directory
cd "$input_dir" || exit 1

# Iterate through files starting with SRR and ending with .fastq.gz
for file in SRR*.fastq; do   
	# Run FastQC on each file and output results to the output directory
	fastqc "$file" -o "$output_dir" --delete
done
