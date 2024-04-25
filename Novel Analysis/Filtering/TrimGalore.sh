# This script performs quality trimming and quality control on paired-end sequencing data using Trim Galore and FastQC.
# It is designed to run as a job on the ALICE high-performance computing environment.
# It sets up job parameters, loads necessary modules, and executes the trimming and QC process across all read pairs in the specified directory.
# Last modified: 25.04.2024

#!/bin/bash

# ALICE Settings
#SBATCH --job-name=trimming_SRP
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=40gb
#SBATCH --time=16:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bg171@student.le.ac.uk
#SBATCH --output=trimming_output_alice

# Load the necessary modules
module load fastqc/0.12.1-hkgpcde

# Define the path to Trim Galore
export PATH=/scratch/alice/b/bg171/TrimGalore-0.6.10:$PATH

# Define the directories for the reads, trimmed output, and FastQC reports
READS_DIR="/scratch/alice/b/bg171/SRP_Data"
TRIMMED_DIR="/scratch/alice/b/bg171/trimmed_data_alice"
TRIMMEDQC_DIR="/scratch/alice/b/bg171/trimmed_fastqc_alice"

# Make sure the output directories exist
mkdir -p ${TRIMMED_DIR}
mkdir -p ${TRIMMEDQC_DIR}

# Loop through each pair of reads and trim them
for file1 in ${READS_DIR}/*_1.fastq
do
    # Construct the file name for the second read in the pair
    file2=${file1/_1.fastq/_2.fastq}

    # Run Trim Galore! and specify the output directories for the trimmed files and FastQC reports
    trim_galore -q 21 --paired --fastqc --fastqc_args "--outdir ${TRIMMEDQC_DIR}" --output_dir ${TRIMMED_DIR} ${file1} ${file2}
done

echo "Trimming complete!"
