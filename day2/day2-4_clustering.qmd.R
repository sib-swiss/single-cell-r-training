# Clustering Script

# Load necessary libraries
library(Seurat)
library(clustree)

# Load dataset
seu <- readRDS("day2/seu_day2-3.rds")

# Find neighbors
seu <- FindNeighbors(seu, dims = 1:25, reduction = "integrated.cca")

# Find clusters with different resolutions
seu <- FindClusters(seu, resolution = seq(0.1, 0.8, by=0.1))

# Display metadata
head(seu@meta.data)

# Visualize cluster tree
clustree(seu@meta.data[,grep("RNA_snn_res", colnames(seu@meta.data))],
         prefix = "RNA_snn_res.")

# Plot UMAP for a selected resolution
DimPlot(seu, group.by = "RNA_snn_res.0.2")

# Exercise: Try different resolutions
DimPlot(seu, group.by = "RNA_snn_res.0.2") | DimPlot(seu, group.by = "RNA_snn_res.0.3") | DimPlot(seu, group.by = "RNA_snn_res.0.4")

# Exercise: Change clustering algorithm
# Available algorithms: 1 (Louvain), 2, 3, 4 (other modularity optimization methods)
seu <- FindClusters(seu, resolution = 0.3)

clustree(seu@meta.data[,grep("RNA_snn_res", colnames(seu@meta.data))],
         prefix = "RNA_snn_res.")


# Save processed dataset
saveRDS(seu, "day2/seu_day2-4.rds")