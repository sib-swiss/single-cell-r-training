#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

#' A scaled, variable-feature-selected fixture, optionally already PCA'd.
make_fixture <- function(n_cells = 300, n_genes = 200, run_pca = TRUE) {
  genes <- sprintf("GENE%04d", seq_len(n_genes))
  cluster <- sample(c("A", "B", "C"), n_cells, replace = TRUE)
  w <- matrix(1, nrow = n_genes, ncol = n_cells, dimnames = list(genes, NULL))
  marker_idx <- split(sample(seq_len(n_genes), 60), rep(c("A", "B", "C"), 20))
  for (cl in names(marker_idx)) {
    w[marker_idx[[cl]], cluster == cl] <- 10
  }
  depth <- round(stats::rlnorm(n_cells, meanlog = log(3000), sdlog = 0.4))
  m <- matrix(0L, nrow = n_genes, ncol = n_cells,
             dimnames = list(genes, sprintf("CELL%04d", seq_len(n_cells))))
  for (j in seq_len(n_cells)) {
    m[, j] <- stats::rmultinom(1, size = depth[j], prob = w[, j] / sum(w[, j]))[, 1]
  }
  seu <- Seurat::CreateSeuratObject(counts = as(m, "CsparseMatrix"), min.cells = 0, min.features = 0)
  seu$true_cluster <- cluster
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 100, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  if (run_pca) seu <- Seurat::RunPCA(seu, npcs = 20, verbose = FALSE)
  seu
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seu <- make_fixture()
  saveRDS(seu, file.path(outdir, "scaled_pca.rds"))
  cat(sprintf("scaled_pca %d genes x %d cells -> scaled_pca.rds\n", nrow(seu), ncol(seu)))
}
