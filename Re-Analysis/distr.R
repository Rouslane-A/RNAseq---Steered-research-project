### Cluster Composition based on cell-cycle states
# load library
library(viridis)

# Extract cell cycle data and cluster labels from the colData slot of your sce object
cell_cycle_data <- sce$cell_cycle 
cluster_labels <- sce$label

# Combine the cluster labels and cell cycle data into a data frame
data_df <- data.frame(cluster = factor(cluster_labels), cell_cycle = cell_cycle_data)

# Calculate the proportion of each cell cycle category within each cluster
cluster_composition <- table(data_df$cluster, data_df$cell_cycle)

# Convert the frequencies to percentages
cluster_composition_percent <- prop.table(cluster_composition, margin = 1) * 100

# Convert the percentages to a data frame
composition_df <- data.frame(cluster = factor(rownames(cluster_composition_percent)),
                             cell_cycle = colnames(cluster_composition_percent),
                             percentage = as.vector(cluster_composition_percent))

# Create custom color palette
my_palette <- inferno(nrow(composition_df))

# Create bar plots using ggplot2
ggplot(composition_df, aes(x = cluster, y = percentage, fill = cell_cycle)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_viridis_d() +  # Apply color scale
  labs(x = "Cluster", y = "Percentage", fill = "Cell Cycle") +
  ggtitle("Composition of Cell Clusters by Cell Cycle") +
  theme_bw()  # Set a white background theme for clarity
