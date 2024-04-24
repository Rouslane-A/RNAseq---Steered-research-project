library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)

counts <- read.csv(file = '/home/bg171/SRP/final_counts.csv', row.names = 1)
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
write.csv(counts, "counts_updated.csv", row.names = TRUE)
