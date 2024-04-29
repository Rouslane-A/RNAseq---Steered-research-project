###initial clustering of all cells based on cell origins
library(scater)
library(scuttle)

#read csv count matrix as dataframe
raw_matrix <- read.csv(file = 'GSE202695_counts_raw.csv', row.names = 1)

#construct SingleCellExperiment cell from count matrix
sce <- SingleCellExperiment(assays = list(counts = raw_matrix))

# Define the pseudocount value
pseudocount <- 1

# Add the pseudocount to the count matrix
counts_with_pseudocount <- counts(sce) + pseudocount

# Replace the original count matrix in the SingleCellExperiment object
assay(sce, "counts") <- counts_with_pseudocount

# Compute log-normalized values for new count matrix
sce <- logNormCounts(sce)

#t-SNE dimensionality reduction using runTSNE
sce <- runTSNE(sce, assay.type = "logcounts")


#Categorize column names based on expressions (cell origins)
categorize_columns <- function(col_name) {
  if (grepl("HBRX3078", col_name)) {
    return("PDX1")
  } else if (grepl("HBRX2353", col_name)) {
    return("PDX2")
  } else if (grepl("HBRX1921", col_name)) {
    return("PDX3")
  } else if (grepl("HBRX2344", col_name)) {
    return("PDX4")
  } else if (grepl("MDAMB231", col_name)) {
    return("MDAMB231")
  } else {
    return("Other")
  }
}


# Apply the function to categorize column names
column_categories <- sapply(colnames(logcounts(sce)), categorize_columns)

# Convert the result to a factor variable
column_categories_factor <- factor(column_categories)

# Add column_categories_factor to colData
colData(sce)$category <- column_categories_factor

# Plot t-SNE plot based on cell origin 
plotTSNE(sce, colour_by = "category")

