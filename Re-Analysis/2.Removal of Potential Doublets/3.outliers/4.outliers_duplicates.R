library(data.table)
library(griph)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Read gene count data from a CSV file
counts_filtered <- fread("/scratch/alice/r/ra500/counts/filtered_counts.csv", header = TRUE)
colnames(counts_filtered)[1] <- "genes"

# Map Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = counts_filtered$genes, keytype = "ENSEMBL", column = 'ENTREZID')
counts_filtered$genes <- entrez_ids

# Convert count data to matrix
count_matrix <- as.matrix(counts_filtered[, -1])

# Set row names as Entrez IDs
rownames(count_matrix) <- counts_filtered$genes

# Remove the prefix from column names
colnames(counts_filtered) <- gsub("_Aligned.out_sorted.bam", "", colnames(counts_filtered))

# Ensure no missing values in the data
clean_counts <- na.omit(counts_filtered)

# Remove duplicate rows
clean_counts <- unique(clean_counts)

# Save filtered_counts as a CSV file
write.csv(clean_counts, file = "/scratch/alice/r/ra500/counts/clean_counts.csv")

# Perform log transformation on the count data
log_counts <- log2(count_matrix + 1)

# Save log_counts as a CSV file
write.csv(log_counts, file = "/scratch/alice/r/ra500/counts/log_counts.csv")

