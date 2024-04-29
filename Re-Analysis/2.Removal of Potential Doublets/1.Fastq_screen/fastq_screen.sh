#!/bin/bash
#
# Example SLURM job script for ALICE

#SBATCH --job-name=features
#SBATCH --cpus-per-task=50
#SBATCH --mem=200gb
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ra500@student.le.ac.uk

module load bowtie2


# Set the path to FastQ Screen executable
FASTQSCREEN_PATH="/home/r/ra500/Downloads/FastQ-Screen-0.15.3/fastq_screen"

# Set the path to the reference database
CONFIG="/home/r/ra500/Downloads/FastQ-Screen-0.15.3/fastq_screen.conf"

# Set the output directory
OUTPUT_DIR="/scratch/alice/r/ra500/fastq_screen_out"

# Function to check if a file has already been processed
is_processed() {
    local file="$1"
    local base_name
    base_name=$(basename "$file" .fastq)
    if [[ -d "$OUTPUT_DIR/$base_name" ]]; then
        return 0  # File has been processed
    else
        return 1  # File has not been processed
    fi
}

# Iterate over each FASTQ file in the folder
for fastq_file in /scratch/alice/r/ra500/rna-seq/*.fastq; do

    # Check if the file has already been processed
    if is_processed "$fastq_file"; then
        echo "Skipping $fastq_file as it has already been processed."
        continue
    fi

    # Get the base name of the file
    base_name=$(basename "$fastq_file" .fastq)
    
    # Run FastQ Screen on the current FASTQ file
    $FASTQSCREEN_PATH --aligner bowtie2 --threads 30 --outdir "$OUTPUT_DIR" --conf "$CONFIG" "$fastq_file"
done

