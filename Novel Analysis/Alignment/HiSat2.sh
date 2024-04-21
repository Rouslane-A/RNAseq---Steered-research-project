#!/bin/bash

# ALICE settings
#SBATCH --job-name=align_SRP
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=40gb
#SBATCH --time=36:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bg171@student.le.ac.uk
#SBATCH --output=align_output

module load hisat2/2.2.1-eovlo7b
module load samtools/1.17-wenuvv5

# Define paths to the appropiate files and directories
FASTQ_DIR="/scratch/alice/b/bg171/trimmed_data"
OUTPUT_DIR="/scratch/alice/b/bg171"
HISAT2_INDEX="/scratch/alice/b/bg171/HiSat2_index/hg38_hisat2_index"
ALIGNMENT_DIR="${OUTPUT_DIR}/aligned_alice"

# Run HISAT2 alignment
for file in ${FASTQ_DIR}/*_1.fastq; do
    base=$(basename "${file}" .1.fastq)
    hisat2 -x "${HISAT2_INDEX}" -1 "${file}" -2 "${FASTQ_DIR}/${base}_2.fastq" -S "${ALIGNMENT_DIR}/${base}.sam"
    samtools view -bS "${ALIGNMENT_DIR}/${base}.sam" | samtools sort -o "${ALIGNMENT_DIR}/${base}_sorted.bam"
    samtools index "${ALIGNMENT_DIR}/${base}_sortedindex.bam"
done

echo "Alignment complete"
