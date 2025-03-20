# Cell annotation

# Load necessary libraries
library(Seurat)
library(celldex)
library(SingleR)
library(dittoSeq)

# Load dataset
seu <- readRDS("day2/seu_day2-4.rds")

# Set identity based on clustering resolution
seu <- Seurat::SetIdent(seu, value = seu$RNA_snn_res.0.3)

# Use original count data
DefaultAssay(seu) <- "RNA"

# Visualize expression of a gene
Seurat::FeaturePlot(seu, "HBA1")

# Define marker genes
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

# Visualize T cell marker expression
Seurat::FeaturePlot(seu, tcell_genes, ncol=2)
Seurat::VlnPlot(seu, features = tcell_genes, ncol = 2)

# Exercise: Check monocyte gene expression
Seurat::FeaturePlot(seu, monocyte_genes, ncol=2)
Seurat::VlnPlot(seu, features = monocyte_genes, ncol = 2)

# Add module score for T cell genes
seu <- Seurat::AddModuleScore(seu, features = list(tcell_genes), name = "tcell_genes")

# Visualize module scores
Seurat::FeaturePlot(seu, "tcell_genes1")
Seurat::VlnPlot(seu, "tcell_genes1")

# Annotate cell cycle phase
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

seu <- Seurat::CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
Seurat::DimPlot(seu, group.by = "Phase")

# Automated annotation using SingleR
ref <- celldex::NovershternHematopoieticData()
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu), ref = ref, labels = ref$label.main)

# Visualize SingleR annotation
SingleR::plotScoreHeatmap(seu_SingleR)
SingleR::plotDeltaDistribution(seu_SingleR)

# Remove low-count annotations
singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

# Add annotations to metadata
seu$SingleR_annot <- singleR_labels

# Plot annotations in UMAP
dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7)
dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "orig.ident")

# Compare manual vs automated annotation
dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "RNA_snn_res.0.3")

# Save dataset
saveRDS(seu, "day3/seu_day3-1.rds")
