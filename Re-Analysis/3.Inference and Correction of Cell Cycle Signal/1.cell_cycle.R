# Load the required libraries
library(data.table)
library(scater)

# Read the log counts data from a CSV file
log_counts <- fread("/scratch/alice/r/ra500/counts/log_counts.csv", header = TRUE)

# Perform cell cycle prediction using the log counts data
cell_cycle <- predictCellCycle(log_counts,
                               org = "human.Whitfield",
                               cor_thr = 0.2,
                               refine_iter = 200)
                               

