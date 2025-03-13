library(Seurat)
library(Matrix)
library(ggplot2)

# vector of paths to all sample directories
datadirs <- list.files(path = "data", full.names = TRUE)

# get the sample names
# replace underscores with hyphen to correctly extract sample names later on
names(datadirs) <- basename(datadirs) |> gsub("_", "-", x = _) 

# for now, we only take the HPV negative and cervical cancer samples
datadirs <- datadirs[c("N-HPV-NEG-1", "N-HPV-NEG-2", "SCC-4", "SCC-5")]

# create a large sparse matrix from all count data
sparse_matrix <- Seurat::Read10X(data.dir = datadirs)

# create a seurat object from sparse matrix
seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "CervicalCancerStudy")

Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))

# mitochondrial genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")

# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")

# hemoglobin genes (but not HBP)
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")



Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))

Seurat::FeatureScatter(seu, 
                       feature1 = "percent.globin", 
                       feature2 = "percent.ribo")


most_expressed_boxplot <- function(object, ngenes = 20){
  
  # matrix of raw counts
  cts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  
  # get percentage/cell
  cts <- t(cts)/colSums(cts)*100
  medians <- SparseArray::colMedians(cts) # apply(cts, 2, median)
  
  # get top n genes
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
  
  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
  
  # boxplot with ggplot2
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}

most_expressed_boxplot(seu, 20)

seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 25)

Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))


# Don't run it yet! Read the exercise first
seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)

seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

seu <- Seurat::ScaleData(seu)

seu <- Seurat::RunPCA(seu)

Seurat::DimPlot(seu, reduction = "pca")

Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.ribo")

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

Seurat::ElbowPlot(seu, ndims = 40)

seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap", group.by = "orig.ident")

Seurat::FeaturePlot(seu, "SIFa")
Seurat::FeaturePlot(seu, "CG10804")



