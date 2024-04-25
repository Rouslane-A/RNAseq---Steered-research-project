# This comprehensive R script utilizes the Seurat package to perform detailed single-cell RNA sequencing analysis.
# Key steps include data importing, preprocessing to ensure compatibility with Seurat, creation of Seurat objects,
# and rigorous quality control checks. It includes normalization, feature selection, and various visualizations such
# as violin plots, scatter plots, and PCA plots to assess data quality and expression patterns. The script proceeds
# with cell cycle scoring, dimensional reduction via PCA and UMAP, and clustering. Additionally, it maps gene symbols
# to Entrez IDs for Gene Set Enrichment Analysis (GSEA), calculates cell cycle phase percentages per cluster,
# identifies and visualizes marker genes, and performs enrichment analysis for clusters. The outputs are formatted
# for insightful biological interpretation and presented through various plots and heatmaps, providing a comprehensive
# toolkit for in-depth analysis of single-cell transcriptomic data.
# Last modified: 25.04.2024

# Importing the packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)

# Import the CSV file
counts <- read.csv('/home/bg171/srp/final_counts.csv', row.names = 1)

## process the data so no rownames or column names have '_' (this is not supported by Seurat)
rownames(counts) <- gsub(pattern = "_", 
                         replacement = "-", 
                         x = rownames(counts))

colnames(counts) <- gsub(pattern = "_", 
                         replacement = "-", 
                         x = colnames(counts))

## create a SeuratObject for analysis
seuratObject <- CreateSeuratObject(counts = counts)

# QUALITY CONTROL
## calculate the percentage of mitochondrial genes
grep("^MT", rownames(seuratObject), value = TRUE)
PercentageFeatureSet(seuratObject, pattern = "^MT") -> seuratObject$percent.MT
head(seuratObject$percent.MT)

## calculate the percentage of ribosomal genes
grep("^RP[LS]",rownames(seuratObject),value = TRUE)
PercentageFeatureSet(seuratObject,pattern="^RP[LS]") -> seuratObject$percent.Ribosomal
head(seuratObject$percent.Ribosomal)

## Visualise the quality control metrics (to help determine filtering thresholds)
VlnPlot(seuratObject, layer="counts", features=c("nCount_RNA","percent.MT"), log = TRUE)

## data visualisation to determine the quality control cut offs
VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"))

##Visualisation of nCount_RNA vs. percent.MT
plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.MT") + 
  theme(legend.position = "none")

## Visualisation of nCount_RNA vs. nFeature_RNA
plot2 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none")
plot1 + plot2

## Filter away cells that have unique feature counts(genes) less than 300 or more than 10000. 
## Also filter away cells that have > 10% mitochondrial counts and more than 3e+06 RNA content
seuratObject <- subset(seuratObject, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & nCount_RNA < 3e+06 & percent.MT < 10)

## Visualising after filtering violin plots
VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3, pt.size = 0.001)

## Visualizing the data after filtering scatter plots
plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.MT") + 
  theme(legend.position = "none")
plot2 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none")
plot1 + plot2

# SEURAT ANALYSIS
## Data Normalisation using log normalisation
seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

## Visualisation of the data normalization
set.seed(123)
par(mfrow=c(1,2))

## original expression distribution
raw_geneExp = GetAssayData(seuratObject, assay = "RNA", slot = "counts") %>% as.vector() %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp, main = "Raw Gene Expression", xlab = "Expression Level", ylab = "Frequency")

## expression distribution after normalization
logNorm_geneExp = GetAssayData(seuratObject, assay = "RNA", slot = "data") %>% as.vector() %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
hist(logNorm_geneExp, main = "Normalized Gene Expression", xlab = "Expression Level", ylab = "Frequency")

## Feature Selection 
seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuratObject), 10)

## Plot variable features with and with labels
plot1 <- VariableFeaturePlot(seuratObject) + 
  theme(legend.position = "top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position = "none")
plot2

## Check which S phase genes are present in your dataset
s.genes <- cc.genes$s.genes[cc.genes$s.genes %in% rownames(seuratObject)]
## Check which G2/M phase genes are present in your dataset
g2m.genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% rownames(seuratObject)]

## score the cell cycle genes
if (length(s.genes) < 10 | length(g2m.genes) < 10) {
  stop("Not enough cell cycle genes present in the dataset to perform scoring.")
}

## Only use the genes that are actually present in the seuratObject
seuratObject <- CellCycleScoring(seuratObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## Visualise the cell cycle scores
VlnPlot(seuratObject, features = c("S.Score", "G2M.Score"), ncol = 2)

## Scaling the data 
all.genes <- rownames(seuratObject)
seuratObject <- ScaleData(seuratObject, features = all.genes, verbose = FALSE)

## Linear dimensional reduction
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject), verbose = FALSE)

## Examine and visualise the PCA results a few different ways
print(seuratObject[["pca"]], dims = 1:5, nfeatures = 5)

## Visuals
VizDimLoadings(seuratObject, dims = 1:2, reduction = "pca")

DimPlot(seuratObject, reduction = "pca")

DimHeatmap(seuratObject, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seuratObject, dims = 1:9, cells = 500, balanced = TRUE)

## Clustering the cells
seuratObject <- FindNeighbors(seuratObject, dims = 1:4, verbose = FALSE)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)

## Look at cluster IDs of the first 5 cells
head(Idents(seuratObject), 5)

## Non-linear dimensional reduction
seuratObject <- RunUMAP(seuratObject, dims = 1:4, verbose = FALSE)

## UMAP plot of the clusters
DimPlot(seuratObject, reduction = "umap")

## tSNE plot of the clusters
seuratObject <- RunTSNE(seuratObject, dims = 1:4, verbose = FALSE)
DimPlot(seuratObject, reduction = "tsne")

## Add cluster labels to the tsne and UMAP plots
DimPlot(seuratObject, reduction = "umap", label = TRUE)

plot <- DimPlot(object = seuratObject)
LabelClusters(plot = plot, id = 'ident')

## Create a new metadata column combining S.Score and G2M.Score
seuratObject$CellCycleStage <- ifelse(seuratObject$S.Score > seuratObject$G2M.Score, 'S', 'G2/M')

## Create a t-SNE and UMAP plot color-coded by the cell cycle stage
DimPlot(seuratObject, reduction = 'tsne', group.by = 'CellCycleStage') + ggtitle("t-SNE by Cell Cycle Stage")
DimPlot(seuratObject, reduction = 'umap', group.by = 'CellCycleStage') + ggtitle("UMAP by Cell Cycle Stage")

## Calculate the percentage of cells in each cell cycle phase for each cluster
cell_cycle_distribution <- seuratObject@meta.data %>%
  group_by(cluster = Idents(seuratObject), CellCycleStage) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

## Create the stacked bar chart using ggplot2 to visualise the cell cycle stages in each cluster
ggplot(cell_cycle_distribution, aes(x = as.factor(cluster), y = Percentage, fill = CellCycleStage)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(x = "Cluster", y = "Percentage") +
  scale_fill_brewer(palette = "Set1") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

## Identify the marker genes for each cluster
markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(head(markers))

# Select top 10 markers for each cluster based on log fold change
top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>%
  ungroup()

## Ensuring the list of genes to pass to DoHeatmap is correct
heatmap_genes <- top10$gene

## Heatmap of top 10 marker genes for each cluster (if any significant markers were found)
heatmap_plot <- DoHeatmap(seuratObject, features = heatmap_genes) + NoLegend()
heatmap_plot <- heatmap_plot + theme(axis.text.y = element_blank()) # Remove row labels for better presentation of the heatmap
print(heatmap_plot)

## DotPlot showing the different expression of genes across the clusters (alternative visualisation to the heatmap)
genes <- top10$gene
DotPlot(seuratObject, features = genes) + 
  theme(  
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 1),  
    axis.text.y = element_text(size = 10)) +
  xlab("Gene") +  
  ylab("Cluster")

#GENE SET ENRICHMENT ANALYSIS (GSEA) 
## Convert to numeric values if they are not already
markers$p_val <- as.numeric(markers$p_val)
markers$avg_log2FC <- as.numeric(markers$avg_log2FC)

## Map gene symbols to Entrez IDs
markers$entrez <- mapIds(org.Hs.eg.db, keys = markers$gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

## Remove NA values that result from the mapping
markers <- markers[!is.na(markers$entrez), ]

## Rank genes based on adjusted p-values (the lower the p-value, the higher the rank)
gene_list_ranks <- setNames(markers$avg_log2FC, markers$entrez)

## Sort the ranks in decreasing order
gene_list_ranks <- sort(gene_list_ranks, decreasing = TRUE)

## Split markers by cluster to create a list of gene vectors, one for each cluster
cluster_markers_list <- split(markers$entrez, markers$cluster)

## Create a function for the GSEA
run_gsea <- function(genes) {
  if (length(genes) == 0) return(NULL)
  pre_ranked_genes <- gene_list_ranks[genes]  
  pre_ranked_genes <- sort(pre_ranked_genes, decreasing = TRUE)  
  gsea_results <- gseKEGG(geneList = pre_ranked_genes, 
                          organism = 'hsa', 
                          nPerm = 30000, 
                          minGSSize = 5, 
                          maxGSSize = 500, 
                          pvalueCutoff = 0.05, 
                          verbose = FALSE)
  return(gsea_results)
}

## Apply the GSEA with the ranked list
gsea_results <- lapply(cluster_markers_list, run_gsea)

## Generate heatmap of the GSEA results
nes_list <- lapply(gsea_results, function(x) {
  if (!is.null(x)) {
    return(as.numeric(x$NES))
  } else {
    return(rep(NA, length(x)))
  }
})

##
nes_matrix <- do.call(cbind, nes_list)
rownames(nes_matrix) <- names(nes_list[[1]])  # Replace this with actual gene set names from your results

## Plot the heatmap
pheatmap(nes_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_legend = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 12,
         display_numbers = FALSE,
         main = "GSEA Heatmap")

## Define the mapping of clusters to superclusters
supercluster_mapping <- c("2" = "A", "4" = "A", "5" = "B", "3" = "B", "7" = "B")

## Create a new column 'superclusters' in the metadata
## Assuming 'supercluster_mapping' is defined appropriately elsewhere in your code
seuratObject@meta.data$superclusters <- supercluster_mapping[as.character(seuratObject@meta.data$seurat_clusters)]

## Define the EMT and proliferation markers
emt_markers <- c("KRT19", "VIM", "SOX9", "FN1", "MYC", "MUC1", "COL1A1", 
                 "MMP2", "KRT14", "ST14", "EPCAM", "TNC", "EMP3", "ZEB1",
                 "ITGB1", "SPARC", "DSG2", "TGFBR1", "KRT18", "MET", "SERPINE1",
                 "ITGAV", "CDH2")
prolif_markers <- c("MCM3", "PCNA", "MKI67")

## Combine markers into one list
specified_genes <- c(emt_markers, prolif_markers)

## Subset Seurat object by superclusters
supercluster_A_cells <- subset(seuratObject, subset = superclusters %in% c("A"))
supercluster_B_cells <- subset(seuratObject, subset = superclusters %in% c("B"))

## Extract the expression data for the specified genes
supercluster_A_expr <- GetAssayData(supercluster_A_cells, assay = "RNA")[specified_genes, ]
supercluster_B_expr <- GetAssayData(supercluster_B_cells, assay = "RNA")[specified_genes, ]

## Ensure the number of columns is equal
num_cells_A <- ncol(supercluster_A_expr)
num_cells_B <- ncol(supercluster_B_expr)
max_num_cells <- max(num_cells_A, num_cells_B)

## Pad with zeros to match the dimensions
if (num_cells_A < max_num_cells) {
  extra_columns_A <- matrix(0, nrow = nrow(supercluster_A_expr), ncol = max_num_cells - num_cells_A)
  supercluster_A_expr <- cbind(supercluster_A_expr, extra_columns_A)
} else if (num_cells_B < max_num_cells) {
  extra_columns_B <- matrix(0, nrow = nrow(supercluster_B_expr), ncol = max_num_cells - num_cells_B)
  supercluster_B_expr <- cbind(supercluster_B_expr, extra_columns_B)
}

## Add a column of NA to create a break between the superclusters
break_column <- matrix(NA, nrow = nrow(supercluster_A_expr), ncol = 1)
colnames(break_column) <- "break"

## Combine the expression matrices for superclusters A and B with a break in between
heatmap_data <- cbind(supercluster_A_expr, break_column, supercluster_B_expr)
rownames(heatmap_data) <- specified_genes
colnames(heatmap_data) <- c(rep("Supercluster A", ncol(supercluster_A_expr)), "Break", rep("Supercluster B", ncol(supercluster_B_expr)))

## Generate heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         breaks = seq(-2, 2, length.out = 101),
         fontsize = 10,
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         main = "Expression of Specified Genes in Superclusters")
