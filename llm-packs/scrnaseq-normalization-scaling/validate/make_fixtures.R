#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

make_fixture <- function(n_cells = 300, n_genes = 300) {
  genes <- sprintf("GENE%04d", seq_len(n_genes))
  depth <- round(stats::rlnorm(n_cells, meanlog = log(3000), sdlog = 0.5))
  m <- matrix(0L, nrow = n_genes, ncol = n_cells,
             dimnames = list(genes, sprintf("CELL%04d", seq_len(n_cells))))
  # A handful of genes get real cell-to-cell variance so there's something
  # for FindVariableFeatures to find; the rest are near-uniform noise.
  var_genes <- sample(genes, 30)
  w <- matrix(1, nrow = n_genes, ncol = n_cells, dimnames = list(genes, NULL))
  cluster <- sample(c("A", "B"), n_cells, replace = TRUE)
  w[var_genes, cluster == "A"] <- 8
  for (j in seq_len(n_cells)) {
    m[, j] <- stats::rmultinom(1, size = depth[j], prob = w[, j] / sum(w[, j]))[, 1]
  }
  Seurat::CreateSeuratObject(counts = as(m, "CsparseMatrix"), min.cells = 0, min.features = 0)
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seu <- make_fixture()
  saveRDS(seu, file.path(outdir, "qc_filtered.rds"))
  cat(sprintf("qc_filtered %d genes x %d cells -> qc_filtered.rds\n", nrow(seu), ncol(seu)))
}
