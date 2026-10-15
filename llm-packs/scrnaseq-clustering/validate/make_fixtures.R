#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

#' A PCA'd fixture with 3 known cell groups, optionally also "integrated"
#' (a second reduction standing in for integrated.cca, plus a planted batch
#' effect in the unintegrated PCA so the two reductions actually differ).
make_fixture <- function(n_cells = 300, n_genes = 200, with_integration = FALSE) {
  genes <- sprintf("GENE%04d", seq_len(n_genes))
  cluster <- sample(c("A", "B", "C"), n_cells, replace = TRUE)
  batch <- sample(c("S1", "S2"), n_cells, replace = TRUE)
  w <- matrix(1, nrow = n_genes, ncol = n_cells, dimnames = list(genes, NULL))
  marker_idx <- split(sample(seq_len(n_genes), 60), rep(c("A", "B", "C"), 20))
  for (cl in names(marker_idx)) w[marker_idx[[cl]], cluster == cl] <- 10
  if (with_integration) {
    batch_genes <- sample(setdiff(seq_len(n_genes), unlist(marker_idx)), 20)
    w[batch_genes, batch == "S1"] <- 8
  }
  depth <- round(stats::rlnorm(n_cells, meanlog = log(3000), sdlog = 0.4))
  m <- matrix(0L, nrow = n_genes, ncol = n_cells,
             dimnames = list(genes, sprintf("CELL%04d", seq_len(n_cells))))
  for (j in seq_len(n_cells)) {
    m[, j] <- stats::rmultinom(1, size = depth[j], prob = w[, j] / sum(w[, j]))[, 1]
  }
  seu <- Seurat::CreateSeuratObject(counts = as(m, "CsparseMatrix"), min.cells = 0, min.features = 0)
  seu$orig.ident <- batch
  seu$true_cluster <- cluster
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 100, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = 15, verbose = FALSE)
  seu <- Seurat::RunUMAP(seu, dims = 1:15, reduction = "pca", verbose = FALSE)
  if (with_integration) {
    # Stand-in "integrated" reduction: same embedding as pca here is enough
    # to exercise the reduction-selection guard; the fixture's job is to
    # test that the pack *chooses* it, not to re-prove integration quality
    # (that's the integration pack's job).
    seu@reductions[["integrated.cca"]] <- seu@reductions[["pca"]]
  }
  seu
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seu <- make_fixture()
  saveRDS(seu, file.path(outdir, "pca_only.rds"))
  cat(sprintf("pca_only %d genes x %d cells -> pca_only.rds\n", nrow(seu), ncol(seu)))
}
