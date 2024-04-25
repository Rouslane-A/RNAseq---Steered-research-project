# This script performs MultiQC on the specifed directories.
# Last modifed: 25.04.2025

#!/bin/bash

# ALICE settings
#SBATCH --job-name=multiqc_SRP
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=40gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bg171@student.le.ac.uk
#SBATCH --output=multiqc_output_alice

# Load neccesary modules
module load python/3.10.12-tggsi7t

# Define the directories where FastQC and fastq_screen results are stored
FASTQC_DIR="/scratch/alice/b/bg171/fastqc_results_alice"

# Define the output directory for MultiQC
MULTIQC_OUTPUT_DIR="/scratch/alice/b/bg171/multiqc_output_alice"

# Create the MultiQC output directory if it doesn't already exist
mkdir -p "${MULTIQC_OUTPUT_DIR}"

# Navigate to the MultiQC output directory
cd "${MULTIQC_OUTPUT_DIR}"

# Run MultiQC on the results
multiqc "${FASTQC_DIR}" -o "${MULTIQC_OUTPUT_DIR}"

echo "MultiQC report generated at ${MULTIQC_OUTPUT_DIR}"
