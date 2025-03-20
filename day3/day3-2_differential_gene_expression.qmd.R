# Differential gene expression

# Load necessary libraries
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("scuttle", quietly = TRUE)) BiocManager::install("scuttle")

library(Seurat)
library(edgeR)
library(limma)
library(dplyr)
library(scuttle)

# Load dataset
seu <- readRDS("day3/seu_day3-1.rds")

# Find all markers for each cluster
de_genes <- Seurat::FindAllMarkers(seu, min.pct = 0.25, only.pos = TRUE)

de_genes <- subset(de_genes, de_genes$p_val_adj < 0.05)
write.csv(de_genes, "day3/de_genes_FindAllMarkers.csv", row.names = FALSE, quote = FALSE)

top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# Generate a dotplot
dittoSeq::dittoDotPlot(seu, vars = unique(top_specific_markers$gene), group.by = "RNA_snn_res.0.3")

# Define T-cell genes
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")

de_genes[de_genes$gene %in% tcell_genes,]

# Differential expression between groups of cells
seu <- Seurat::SetIdent(seu, value = "SingleR_annot")

deg_cd8_cd4 <- Seurat::FindMarkers(seu,
                                   ident.1 = "CD8+ T cells",
                                   ident.2 = "CD4+ T cells",
                                   group.by = seu$SingleR_annot,
                                   test.use = "wilcox")

deg_cd8_cd4 <- subset(deg_cd8_cd4, deg_cd8_cd4$p_val_adj < 0.05)

# Check CD8A, CD8B, and CD4 expression
deg_cd8_cd4[c("CD4", "CD8A", "CD8B"),]

# Violin plot of genes
Seurat::VlnPlot(seu, features = c("CD4", "CD8A", "CD8B"), idents = c("CD8+ T cells", "CD4+ T cells"))

# Differential expression using limma
proB <- readRDS("course_data/proB.rds")

Seurat::DimPlot(proB, group.by = "orig.ident")

table(proB@meta.data$type)
head(proB@meta.data)

# Prepare pseudobulk count matrix
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

proB$patient.id <- gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id <- sapply(strsplit(proB$patient.id, "-"), '[', 2)

proB$sample <- factor(proB$orig.ident)

bulk <- Seurat::AggregateExpression(proB, group.by = "sample", return.seurat = TRUE, assay = "RNA")

meta_data <- unique(proB@meta.data[, c("orig.ident", "sample", "type", "patient.id")])
rownames(meta_data) <- meta_data$orig.ident
bulk@meta.data <- meta_data[colnames(bulk), ]

counts <- Seurat::GetAssayData(bulk, layer = "counts") |> as.matrix()

y <- edgeR::DGEList(counts, samples = bulk@meta.data)
keep <- edgeR::filterByExpr(y, group = bulk$type)
y <- y[keep,]

# Generate design matrix
design <- model.matrix(~0 + y$samples$type + y$samples$patient.id)
colnames(design) <- make.names(c("ETV6-RUNX1", "PBMMC", "patient2", "patient3"))
rownames(design) <- rownames(y$samples)

# Define contrast
contrast.mat <- limma::makeContrasts(ETV6.RUNX1 - PBMMC, levels = design)

dge <- edgeR::calcNormFactors(y)
vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)

# Get top differentially expressed genes
limma::topTable(fit.contrasts, number = 10, sort.by = "P")
limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val < 0.05))

# Violin plot for selected genes
Seurat::VlnPlot(proB, "S100A9", split.by = "type")
Seurat::VlnPlot(proB, "SOCS2", split.by = "type")

# Run differential analysis with Seurat (without paired design)
tum_vs_norm <- Seurat::FindMarkers(proB, 
                                   ident.1 = "ETV6-RUNX1", 
                                   ident.2 = "PBMMC", 
                                   group.by = "type")
tum_vs_norm <- subset(tum_vs_norm, tum_vs_norm$p_val_adj < 0.05)

# Merge results and compare
merge_limma_FindMarkers <- merge(tum_vs_norm, limma_de, by = "row.names", all.x = TRUE)

par(mar = c(4, 4, 4, 4))
plot(merge_limma_FindMarkers$avg_log2FC, merge_limma_FindMarkers$logFC,
     xlab = "log2FC Wilcoxon", ylab = "log2FC limma",
     pch = 15, cex = 0.5)
abline(a = 0, b = 1, col = "red")

# Save results
saveRDS(tum_vs_norm, "day3/tum_vs_norm_day3-2.rds")
saveRDS(limma_de, "day3/limma_de_day3-2.rds")
saveRDS(proB, "day3/proB_day3-2.rds")