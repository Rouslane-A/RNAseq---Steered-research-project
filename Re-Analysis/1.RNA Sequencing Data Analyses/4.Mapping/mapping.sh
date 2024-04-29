#!/bin/bash

# Alice settings

#SBATCH --job-name=STAR_SRA
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=200gb
#SBATCH --time=23:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ra500@student.le.ac.uk
#SBATCH --output=output_STAR-smp_4


# Define directories
REFERENCE_DIR="/scratch/alice/r/ra500/hg38"
INDEX_DIR="/scratch/alice/r/ra500/seq"
ALIGN_DIR="/scratch/alice/r/ra500/alignment"
FASTQ_DIR="/scratch/alice/r/ra500/rna-seq"
PROCESSED_FILE="$ALIGN_DIR/processed_files.txt"

# Load STAR module
module load star/2.7.10b-m3zkpic


# Check if the processed files list exists
if [ ! -f "$PROCESSED_FILE" ]; then
    touch "$PROCESSED_FILE"
fi

# Perform alignment for each pair of FASTQ files
for fastq1 in "$FASTQ_DIR"/*_1.fastq; do
    fastq2="${fastq1%_1.fastq}_2.fastq"
    sample_id=$(basename "$fastq1" _1.fastq)
    
    # Check if the file has been processed
    if grep -q "^$sample_id$" "$PROCESSED_FILE"; then
        echo "Skipping already processed file: $fastq1"
        continue
    fi

    # Run STAR alignment
    STAR --runMode alignReads \
         --runThreadN 30 \
         --genomeDir "$INDEX_DIR" \
         --readFilesIn "$fastq1" "$fastq2" \
         --sjdbGTFfile "$REFERENCE_DIR/hg38.refGene.gtf" \
         --sjdbOverhang 100 \
         --outFileNamePrefix "$ALIGN_DIR/${sample_id}_" \
         --outSAMstrandField intronMotif \
         --outFilterIntronMotifs RemoveNoncanonical

    # Add the processed file to the list
    echo "$sample_id" >> "$PROCESSED_FILE"
done

