# Load required library
library(Seurat)

# Load the dataset
seu <- readRDS("day1/seu_day1-2.rds")

### Normalization

# Before normalization: inspect raw counts
Seurat::GetAssayData(seu)[1:10,1:10]  

# After normalization
seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)

Seurat::GetAssayData(object = seu, slot = "data")[1:10,1:10]  

# Updating `seu`
# The output of the normalization function is added to the object without losing information

### Variable features

# Identify highly variable features
seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

# Plot the top variable features
vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

# Ensure the plotting window is large enough to avoid errors

### Scaling

# Scale the data
seu <- Seurat::ScaleData(seu)

### The use of `Seurat::SCTransform`

# The functions NormalizeData, VariableFeatures, and ScaleData can be replaced with SCTransform.
# SCTransform performs normalization and scaling in a more sophisticated manner but is slower.

# Bonus exercise: Run SCTransform
seu <- Seurat::SCTransform(seu)

# Check where the output is stored
names(seu@assays)

# Running SCTransform will change `@active.assay` to `SCT` instead of `RNA`
# Check the active assay
DefaultAssay(seu)

# To change the active assay back to `RNA`
DefaultAssay(seu) <- "RNA"

### Save the dataset and clear environment

# Save the dataset
saveRDS(seu, "day2/seu_day2-1.rds")


