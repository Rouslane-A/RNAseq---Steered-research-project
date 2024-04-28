# This R script performs multiple data manipulation steps on RNA sequencing counts data.
# It converts Ensembl gene IDs to gene symbols, maps SRR identifiers to GSM codes,
# and further maps GSM codes to sample titles for enhanced readability and downstream analysis.
# The script utilizes packages such as AnnotationDbi, dplyr, readr, and stringr to efficiently handle and transform the data.
# Outputs are systematically written to updated CSV files, reflecting each stage of transformation within a biological data analysis workflow.
# Last modified: 25.04.2024

library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
library(dplyr)
library(readr)
library(stringr)
library(ensembldb)

# CONVERTING THE ENSEMBLE IDS TO GENE SYMBOLS
counts <- read.csv(file = '/home/bg171/SRP/counts_backup.csv', row.names = 1)
ensembl_ids <- rownames(counts)
head(ensembl_ids)

# Retrieve gene symbols using ensembldb
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl_ids, 
                              keytype = 'GENEID', columns = c("SYMBOL"))

# Create a named vector of gene symbols with Ensembl IDs as names 
gene_symbols <- setNames(geneIDs1$SYMBOL, geneIDs1$GENEID)

# Match the gene symbols with the Ensembl IDs in the counts matrix row names
matched_gene_symbols <- gene_symbols[rownames(counts)]

# Create unique identifiers for duplicated gene symbols by appending Ensembl IDs
duplicated_symbols <- duplicated(matched_gene_symbols) | duplicated(matched_gene_symbols, fromLast = TRUE)
matched_gene_symbols[duplicated_symbols] <- paste(matched_gene_symbols[duplicated_symbols], rownames(counts)[duplicated_symbols], sep = "_")

# Replace the row names in the counts matrix with the unique gene symbols
rownames(counts) <- matched_gene_symbols
head(counts)

# Write the updated counts dataframe to a new file
write.csv(counts, "/home/bg171/srp/counts_updated.csv", row.names = TRUE)

# CONVERTING THE SRR FILE NAMES TO GSM CODES
# Load the SraRunTable to map SRR to GSM
sra_run_table <- read_csv("/home/bg171/SRP/SraRunTable.txt")
srr_to_gsm <- setNames(sra_run_table$`Sample Name`, sra_run_table$Run)

# Load the counts_updated.csv
counts_updated <- read_csv("/home/bg171/srp/counts_updated.csv")

# Update the columns from SRR to GSM
updated_columns <- sapply(strsplit(colnames(counts_updated), "_", fixed = TRUE), function(x) srr_to_gsm[x[1]] %>% coalesce(x[1]))
colnames(counts_updated) <- updated_columns

# Save the updated DataFrame with GSM codes
write_csv(counts_updated, "/home/bg171/srp/counts_updated_with_GSMs.csv")

print("SRR to GSM conversion complete. The updated file is saved as 'counts_updated_with_GSMs.csv'.")

# CONVERTING THE GSM CODES TO THE SAMPLE TITLES
# Read the series matrix file to create the mapping
file_path <- '/home/bg171/SRP/GSE202695_series_matrix.txt'
file <- readLines(file_path)

# Extract GSM codes and sample titles
gsm_codes <- strsplit(file[35], "\t")[[1]]
sample_titles <- strsplit(file[34], "\t")[[1]]

# Remove the first element which is the identifier
gsm_codes <- gsm_codes[-1]
sample_titles <- sample_titles[-1]

# Clean up GSM codes and sample titles by removing quotes
gsm_codes <- gsub('"', '', gsm_codes)
sample_titles <- gsub('"', '', sample_titles)

# Create the mapping from GSM codes to sample titles
gsm_to_sample_title <- setNames(sample_titles, gsm_codes)

# Load the counts dataframe
counts_file_path <- '/home/bg171/srp/counts_updated_with_GSMs.csv'
counts_updated_with_gsms <- read_csv(counts_file_path)

# Update the column names in the counts dataframe using the mapping
colnames(counts_updated_with_gsms) <- sapply(colnames(counts_updated_with_gsms), function(gsm) {
  if (gsm %in% names(gsm_to_sample_title)) {
    return(gsm_to_sample_title[gsm])
  } else {
    return(gsm)  # Retain the original GSM code if it's not found in the mapping
  }
})

# Save the updated dataframe to a new CSV file
write_csv(counts_updated_with_gsms, '/home/bg171/srp/final_counts.csv')
