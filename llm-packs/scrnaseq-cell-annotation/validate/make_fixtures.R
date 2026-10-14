#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

#' A clustered query object with 3 true cell types, plus a small synthetic
#' SingleR reference with 4 profiles (3 matching + 1 "decoy" type absent from
#' the query, so the "every cell still gets a label" behaviour is testable
#' without a real celldex download).
make_fixture <- function(n_cells = 240, n_genes = 120) {
  genes <- sprintf("GENE%04d", seq_len(n_genes))
  types <- c("Tcell", "Bcell", "Monocyte")
  cluster <- sample(types, n_cells, replace = TRUE)
  marker_idx <- split(sample(seq_len(n_genes), 60), rep(types, 20))

  w <- matrix(1, nrow = n_genes, ncol = n_cells, dimnames = list(genes, NULL))
  for (ty in types) w[marker_idx[[ty]], cluster == ty] <- 10
  depth <- round(stats::rlnorm(n_cells, meanlog = log(3000), sdlog = 0.3))
  m <- matrix(0L, nrow = n_genes, ncol = n_cells,
             dimnames = list(genes, sprintf("CELL%04d", seq_len(n_cells))))
  for (j in seq_len(n_cells)) {
    m[, j] <- stats::rmultinom(1, size = depth[j], prob = w[, j] / sum(w[, j]))[, 1]
  }
  seu <- Seurat::CreateSeuratObject(counts = as(m, "CsparseMatrix"), min.cells = 0, min.features = 0)
  seu$true_type <- cluster
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)

  # Reference: mean profile per type, plus a decoy type with its own marker
  # genes the query was never given any cells of.
  decoy_genes <- sample(setdiff(seq_len(n_genes), unlist(marker_idx)), 20)
  ref_types <- c(types, "NK_cell_DECOY")
  ref_mat <- matrix(1, nrow = n_genes, ncol = length(ref_types),
                    dimnames = list(genes, ref_types))
  for (ty in types) ref_mat[marker_idx[[ty]], ty] <- 50
  ref_mat[decoy_genes, "NK_cell_DECOY"] <- 50
  ref_mat <- log1p(ref_mat)

  list(seu = seu, ref = ref_mat, labels = colnames(ref_mat))
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fx <- make_fixture()
  saveRDS(fx$seu, file.path(outdir, "clustered.rds"))
  saveRDS(list(ref = fx$ref, labels = fx$labels), file.path(outdir, "reference.rds"))
  cat(sprintf("clustered %d genes x %d cells -> clustered.rds; reference with %d types -> reference.rds\n",
              nrow(fx$seu), ncol(fx$seu), length(fx$labels)))
}
