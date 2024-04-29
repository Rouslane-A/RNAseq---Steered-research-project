library(Rsubread)

# Set your working directory where your BAM files are located
setwd("/scratch/alice/r/ra500/sorted_indexed/")

# Define paths to your BAM files
bam_files <- list.files(pattern = ".bam$", full.names = TRUE)

# Specify the annotation file (GTF or GFF format) for feature counting
annotation_file <- "/scratch/alice/r/ra500/hg38/Homo_sapiens.GRCh38.111.gtf"

# Define output directory for feature counts
output_dir <- "/scratch/alice/r/ra500/hg38/"

# Create output directory if it doesn't exist
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

# Initialize an empty data frame to store combined counts
combined_counts <- NULL

# Loop through each BAM file and perform feature counting
for (bam_file in bam_files) {
  # Generate output file path based on input BAM file name
  output_file <- paste0(output_dir, "combined_feature_counts.txt")
  
  # Run featureCounts
  counts <- featureCounts(files = bam_file,
                          annot.ext = annotation_file,
                          isGTFAnnotationFile = TRUE,
                          isPairedEnd = TRUE, # Change to TRUE if your data is paired-end
                          countMultiMappingReads = TRUE,
                          nthreads = 30, # Number of threads to use
                          verbose = TRUE)
  
  # Extract count data
  count_data <- counts$counts
  
  # Add a column with sample ID (file name without extension)
  sample_id <- gsub(".bam", "", basename(bam_file))
  count_data$SampleID <- sample_id
  
  # If combined_counts is NULL, set it to count_data
  if (is.null(combined_counts)) {
    combined_counts <- count_data
  } else {
    # Otherwise, merge count_data with combined_counts
    combined_counts <- rbind(combined_counts, count_data)
  }
  
  # Write combined counts to output file after processing each BAM file
  write.table(combined_counts, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("Feature counting completed for", sample_id, "\n")
}
