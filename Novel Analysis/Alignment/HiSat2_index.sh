#!/bin/bash

# Load the HiSat2 module
module load hisat2

# Build the HISAT2 index
echo "Building the HISAT2 index for hg38..."
hisat2-build hg38.fa hg38_hisat2_index

echo "Indexing complete. HISAT2 index files for hg38 have been created."
