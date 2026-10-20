#!/usr/bin/env Rscript
#
# Synthetic fixtures for CI. Deliberately NOT derived from course_data - fast,
# deterministic, tiny, and shaped specifically to make the two failure modes
# this pack exists to catch reproducible:
#
#   1. A technical batch effect that integration SHOULD remove.
#   2. A biological ("condition") effect, NOT confounded with batch, that
#      integration should NOT remove - this is what
#      check_integration_preserves_signal() in implementation/integrate.R is
#      tested against.
#
# Realism is only needed in the rough shape of the count distribution and in
# there being a genuine, detectable structure to find or lose - not in the
# biology.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

#' A two-batch, two-condition fixture, condition crossed with batch (NOT
#' confounded): each batch contains cells from both conditions.
make_integration_fixture <- function(n_cells_per_group = 75, n_other = 150) {
  n_batch_genes <- 20
  n_cond_genes  <- 20
  batch_genes <- sprintf("BATCHGENE%02d", seq_len(n_batch_genes))
  cond_genes  <- sprintf("CONDGENE%02d", seq_len(n_cond_genes))
  other_genes <- sprintf("GENE%03d", seq_len(n_other))
  genes <- c(batch_genes, cond_genes, other_genes)

  groups <- expand.grid(batch = c("S1", "S2"), condition = c("healthy", "disease"),
                        stringsAsFactors = FALSE)
  cells <- list(); meta <- list()
  for (g in seq_len(nrow(groups))) {
    batch <- groups$batch[g]; cond <- groups$condition[g]
    n <- n_cells_per_group
    w <- matrix(1, nrow = length(genes), ncol = n,
                dimnames = list(genes, NULL))
    # Technical batch effect: shifts batch_genes, independent of condition.
    w[batch_genes, ] <- if (batch == "S1") 8 else 1
    # Biological condition effect: shifts cond_genes, independent of batch.
    w[cond_genes, ]  <- if (cond == "disease") 8 else 1

    depth <- round(rlnorm(n, meanlog = log(2000), sdlog = 0.4))
    m <- matrix(0L, nrow = length(genes), ncol = n, dimnames = list(genes, NULL))
    for (j in seq_len(n)) {
      m[, j] <- rmultinom(1, size = depth[j], prob = w[, j] / sum(w[, j]))[, 1]
    }
    colnames(m) <- sprintf("%s_%s_CELL%03d", batch, cond, seq_len(n))
    cells[[g]] <- m
    meta[[g]] <- data.frame(orig.ident = rep(batch, n), condition = rep(cond, n),
                            row.names = colnames(m), stringsAsFactors = FALSE)
  }

  counts <- do.call(cbind, cells)
  md <- do.call(rbind, meta)

  seu <- Seurat::CreateSeuratObject(counts = as(counts, "CsparseMatrix"),
                                    meta.data = md, min.cells = 0, min.features = 0)
  seu
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seu <- make_integration_fixture()
  saveRDS(seu, file.path(outdir, "integration.rds"))
  cat(sprintf("integration %d genes x %d cells (2 batches x 2 conditions) -> integration.rds\n",
              nrow(seu), ncol(seu)))
}
