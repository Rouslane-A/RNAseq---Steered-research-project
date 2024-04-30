### Graph-Based Clustering & Gene-set enrichment analyses

# Install packages for GSEA
install.packages("msigdbr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
# Relevant libraries for GSEA
library(scran)
library(scater)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(BiocGenerics)

# Perform graph-based clustering on log-normalized counts
g <- buildSNNGraph(sce)
cluster <- igraph::cluster_walktrap(g)$membership

# Assigning to the 'colLabels' of the SCE object
colLabels(sce) <- factor(cluster)
table(colLabels(sce))
sce <- runTSNE(sce)
plotTSNE(sce, colour_by="label")

# Find markers based on graph-based clustering
out <- findMarkers(sce, pval.type="some", min.prop=0.2, groups = sce$label, 
                   test.type = "t")

# Selected clusters for plotting
cluster1 <- out[[1]]
cluster8 <- out[[8]]
cluster24 <- out[[24]]

# Calculate normalized marker effect sizes
logFC1 <- getMarkerEffects(cluster1)
logFC8 <- getMarkerEffects(cluster8)
logFC24 <- getMarkerEffects(cluster24)

#subset markers to get top 25 genes with largest normalized effect sizes and plot as heatmaps
pheatmap(logFC1[1:25, ], breaks=seq(-5, 5, length.out=101))
pheatmap(logFC8[1:25, ], breaks=seq(-5, 5, length.out=101))
pheatmap(logFC24[1:25, ], breaks=seq(-5, 5, length.out=101))
