#!/bin/bash

# ALICE settings
#SBATCH --job-name=counts
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150gb
#SBATCH --time=18:00:00
#SBATCH --output=$(pwd)/counts_%j.out
#SBATCH --error=$(pwd)/counts_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bg171@student.le.ac.uk
#SBATCH --export=NONE

# Get the current working directory
work_dir=$(pwd)
module load R/4.3.1
R --vanilla << EOF

# Check and install required packages if not already installed
if (!requireNamespace("ensembldb", quietly = TRUE)) {
 if (!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager")
 }
 BiocManager::install("ensembldb")
}
if (!requireNamespace("Rsubread", quietly = TRUE)) {
 if (!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager")
 }
 BiocManager::install("Rsubread")
}

# Load required packages
library(Rsubread)
library(ensembldb)

# Set input and output paths relative to the current working directory
bam_dir <- file.path("$work_dir", "aligned_alice")
output_dir <- file.path("$work_dir", "counts")
ref_dir <- file.path("$work_dir", "annotation")

# Get sorted BAM files
bam_files <- list.files(path=bam_dir, pattern="*_sorted.bam$", full.names=TRUE)

# Download and prepare the GTF file
gtf_file <- file.path(ref_dir, "Homo_sapiens.GRCh38.111.gtf")

# Download the GTF file and place it in the ref_dir
# Run featureCounts
bamcounts <- featureCounts(
 files = bam_files,
 annot.ext = gtf_file,
 isGTFAnnotationFile = TRUE,
 GTF.featureType = "exon",
 GTF.attrType = "gene_id",
 useMetaFeatures = TRUE,
 allowMultiOverlap = TRUE,
 isPairedEnd = TRUE,
 requireBothEndsMapped = TRUE,
 countChimericFragments = FALSE,
 nthreads = 30
)

# print structure and dimensions of bamcounts and its components
print(str(bamcounts))
print(dim(bamcounts$counts))
print(dim(bamcounts$annotation))
print(dim(bamcounts$targets))
print(dim(bamcounts$stat))

# Write counts matrix to file
for (n in names(bamcounts)){
  write.table(bamcounts[[n]], file=file.path(output_dir, paste0(n, ".csv")), sep=",", quote=F, col.names=NA)
}

EOF

