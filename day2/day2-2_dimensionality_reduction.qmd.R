# Load required package
library(Seurat)

# Load dataset
seu <- readRDS("day2/seu_day2-1.rds")

# Run PCA
seu <- Seurat::RunPCA(seu)

# View PCA plot
Seurat::DimPlot(seu, reduction = "pca", dims = c(13,44))

# Color PCA plot by a feature
Seurat::FeaturePlot(seu, reduction = "pca", 
                    features = "JCHAIN", dims = c(2,3))

# Identify top variable genes
head(VariableFeatures(seu))

# Generate heatmap for top principal components
Seurat::DimHeatmap(seu, dims = 1:4, cells = 50, balanced = TRUE, nfeatures = 5)

# Generate elbow plot
Seurat::ElbowPlot(seu, ndims = 40) +
  ggplot2::geom_vline(xintercept = 25, linetype = "dashed", color = "red", linewidth = 1)

# Run UMAP with default settings\
seu <- Seurat::RunUMAP(seu, dims = 1:25)

# View UMAP plot
Seurat::DimPlot(seu, reduction = "umap")

# Color UMAP plot by features
Seurat::FeaturePlot(seu, reduction = "pca", features = c("HBA1", "DNM1L", "percent.globin", "IGKC", "percent.mito"))

# Change UMAP parameters
seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
Seurat::DimPlot(seu, reduction = "umap")

# Test with different number of PCs
seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = "umap")

set.seed(35)
seu <- Seurat::RunUMAP(seu, dims = 1:50)
c50 <- Seurat::DimPlot(seu, reduction = "umap")

set.seed(35)
# Restore UMAP with default PCs for consistency
seu <- Seurat::RunUMAP(seu, dims = 1:25)
c25 <- Seurat::DimPlot(seu, reduction = "umap")

library(patchwork)

(c25 | c25) / c50

# Save dataset
saveRDS(seu, "day2/seu_day2-2.rds")
