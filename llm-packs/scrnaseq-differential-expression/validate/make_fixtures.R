#!/usr/bin/env Rscript
#
# Synthetic fixtures for CI. Deliberately NOT derived from course_data - fast,
# deterministic, tiny, and shaped specifically to make the pseudoreplication
# failure mode (rule R-DE-001) reproducible in a controlled way:
#
#   - TRUEDE genes: a consistent shift between conditions, present in every
#     patient - both a naive per-cell test and pseudobulk should find these.
#   - NOISE genes: a strong patient-to-patient baseline shift that is
#     UNRELATED to condition (pure "some patients just run hotter for this
#     gene" variation) - a naive per-cell test, with hundreds of cells per
#     patient amplifying false precision, is expected to call many of these
#     "significant" by chance; pseudobulk, which has only as many data points
#     as there are patients, should not.
#
# This lets validate/run_all.R reproduce the course's own
# 3820-genes-vs-far-fewer contrast on a tiny, fast, controlled dataset.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

#' @param conditions         character vector of condition labels (2 for the
#'                           standard pseudobulk_de() test, 3+ for anova_de()).
#' @param n_samples_per_cond independent replicates (patients) per condition.
#' @param n_cells_per_sample cells contributed by each patient.
make_de_fixture <- function(conditions = c("A", "B"), n_samples_per_cond = 3,
                            n_cells_per_sample = 60, n_other = 150) {
  n_truede <- 10; n_noise <- 10
  truede_genes <- sprintf("TRUEDE%02d", seq_len(n_truede))
  noise_genes  <- sprintf("NOISE%02d",  seq_len(n_noise))
  other_genes  <- sprintf("GENE%03d",  seq_len(n_other))
  genes <- c(truede_genes, noise_genes, other_genes)

  samples <- expand.grid(condition = conditions, rep = seq_len(n_samples_per_cond),
                         stringsAsFactors = FALSE)
  samples$sample_id <- sprintf("%s-P%d", samples$condition, samples$rep)
  # Every patient gets its own random noise-gene multiplier, independent of
  # condition - this is the "patient effect" that pseudoreplication mistakes
  # for a condition effect when pooled naively at the cell level.
  samples$noise_mult <- stats::runif(nrow(samples), 0.3, 6)

  cells <- list(); meta <- list()
  for (g in seq_len(nrow(samples))) {
    cond <- samples$condition[g]; sid <- samples$sample_id[g]
    n <- n_cells_per_sample
    w <- matrix(1, nrow = length(genes), ncol = n, dimnames = list(genes, NULL))
    # True, condition-consistent effect: present for every patient of that
    # condition (small per-cell jitter only).
    w[truede_genes, ] <- (match(cond, conditions) * 4) * stats::runif(n, 0.85, 1.15)
    # Patient-specific noise, unrelated to condition.
    w[noise_genes, ] <- samples$noise_mult[g] * stats::runif(n, 0.85, 1.15)

    depth <- round(stats::rlnorm(n, meanlog = log(1500), sdlog = 0.3))
    m <- matrix(0L, nrow = length(genes), ncol = n, dimnames = list(genes, NULL))
    for (j in seq_len(n)) {
      m[, j] <- stats::rmultinom(1, size = depth[j], prob = w[, j] / sum(w[, j]))[, 1]
    }
    colnames(m) <- sprintf("%s_CELL%03d", sid, seq_len(n))
    cells[[g]] <- m
    meta[[g]] <- data.frame(orig.ident = rep(sid, n), patient.id = rep(sid, n),
                            type = rep(cond, n),
                            row.names = colnames(m), stringsAsFactors = FALSE)
  }

  counts <- do.call(cbind, cells)
  md <- do.call(rbind, meta)
  Seurat::CreateSeuratObject(counts = as(counts, "CsparseMatrix"), meta.data = md,
                            min.cells = 0, min.features = 0)
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seu2 <- make_de_fixture(c("A", "B"))
  saveRDS(seu2, file.path(outdir, "two_condition.rds"))
  cat(sprintf("two_condition %d genes x %d cells -> two_condition.rds\n", nrow(seu2), ncol(seu2)))
  seu3 <- make_de_fixture(c("young", "middle", "old"))
  saveRDS(seu3, file.path(outdir, "three_condition.rds"))
  cat(sprintf("three_condition %d genes x %d cells -> three_condition.rds\n", nrow(seu3), ncol(seu3)))
}
