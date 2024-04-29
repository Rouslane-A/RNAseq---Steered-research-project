# Load the required libraries
library(scater)
library(griph)
library("org.Hs.eg.db")
library("AnnotationDbi")

# Extract log-normalized counts from the SCE object
log_matrix <- logcounts(sce)

# Extract corresponding Entrez ID names for Gene symbols in counts
entrez_symbol <- select(org.Hs.eg.db, keys = rownames(log_matrix), column = "ENTREZID", 
                        keytype = "SYMBOL")


# Drop values where an Entrez ID could not be found for a Gene symbol
entrez_symbol_clean <- na.omit(entrez_symbol)

# Using cleaned values; create dictionary of Entrez IDs and their Gene symbols
gene_symbol_to_entrez <- setNames(entrez_symbol_clean$ENTREZID, entrez_symbol_clean$SYMBOL)

# Subset log matrix based on gene symbols in the dictionary
sublog_matix <- log_matrix[rownames(log_matrix) %in% names(gene_symbol_to_entrez), ]

# Replace gene symbols in the row names of sublog_matrix with corresponding Entrez IDs from the dictionary
rownames(sublog_matix) <- gene_symbol_to_entrez[rownames(sublog_matix)]

# Perform cell cycle prediction using the log counts data
cell_cycle <- predictCellCycle(sublog_matix,
                               org = "human.Whitfield",
                               cor_thr = 0.2,
                               refine_iter = 200)


# Create new dataframe with cell names and corresponding cell cycle state
cc_dataframe <- data.frame(row.names = TRUE, colnames(sce), cell_cycle)

# Assign dataframe into colData slot of the SCE object
colData(sce) <- cbind(colData(sce), cc_dataframe)

# t-SNE dimensionality reduction using runTSNE
sce <- runTSNE(sce, assay.type = "logcounts")

# Plot t-SNE plot based on cell origin 
plotTSNE(sce, colour_by = "cell_cycle")


