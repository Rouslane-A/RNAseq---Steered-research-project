#!/bin/bash

# Alice settings

#SBATCH --job-name=STAR_SRA
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=40gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ra500@student.le.ac.uk
#SBATCH --output=output_STAR-smp_4

# loading the module sratoolkit
module load sratoolkit/3.0.0-5fetwpi

# Path to the file containing the list of SRR accessions downloaded from the SRA website from the GEO
accession_list="SRR_Acc_List.txt"

# Loop through each accession in the list, 
while IFS= read accession; do
	prefetch $accession  # Download SRA data
	fasterq-dump --progress --split-files $accession  # Convert to FASTQ
    	rm -rf $accession
done < "$accession_list"


