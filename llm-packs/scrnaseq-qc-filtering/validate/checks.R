#!/usr/bin/env Rscript
#
# Runtime post-conditions from rules/qc-filtering.rules.yml. Needs Seurat, so
# this runs in the course container, not on a bare R image.
#
# Two uses:
#   1. CI, against validate/fixtures/*.rds (see run_all.R)
#   2. An agent checking its own output: source this, call check_qc_object(seu)
#      and report every check it could not confirm.

suppressPackageStartupMessages(library(Seurat))

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on an object that has had QC metrics added
check_qc_object <- function(seu) {
  res <- list()

  # R-001: percent.mito identically zero is a pattern/species mismatch.
  pm <- seu$percent.mito
  res[[length(res)+1]] <- chk(
    "C-001", "R-001",
    !is.null(pm) && max(pm, na.rm = TRUE) > 0,
    "percent.mito must not be 0 for every cell (indicates an empty pattern match)")

  # R-001b: an implausibly high median suggests the pattern is over-broad.
  res[[length(res)+1]] <- chk(
    "C-002", "R-001",
    !is.null(pm) && stats::median(pm, na.rm = TRUE) < 50,
    "median percent.mito >= 50% suggests an over-broad mitochondrial pattern")

  # R-004: counts must be joined before any single-matrix operation.
  layers <- tryCatch(SeuratObject::Layers(seu[["RNA"]]), error = function(e) character(0))
  n_count_layers <- sum(grepl("^counts", layers))
  res[[length(res)+1]] <- chk(
    "C-004", "R-004", n_count_layers <= 1,
    sprintf("%d counts layers present; call JoinLayers() before operating on a single counts matrix",
            n_count_layers))

  # R-013: QC must precede integration.
  red <- names(seu@reductions)
  res[[length(res)+1]] <- chk(
    "C-013", "R-013",
    !length(red) || !any(grepl("integrated|harmony|cca|mnn", red, ignore.case = TRUE)),
    "object carries an integrated reduction; QC filtering belongs before integration")

  do.call(rbind, res)
}

#' Post-conditions on the result of filter_cells()
check_filtered_object <- function(seu_before, seu_after) {
  res <- list()
  frac <- 1 - ncol(seu_after) / ncol(seu_before)

  # R-011
  res[[length(res)+1]] <- chk(
    "C-011", "R-011", frac < 0.40,
    sprintf("removed %.1f%% of cells (guard: 40%%)", 100 * frac))

  res[[length(res)+1]] <- chk(
    "C-014", "R-011", ncol(seu_after) > 0,
    "filtering must not remove every cell")

  # The filter must actually have been applied.
  res[[length(res)+1]] <- chk(
    "C-015", "R-010", frac > 0,
    "no cells removed - verify the thresholds were applied and are meaningful")

  do.call(rbind, res)
}

report <- function(res, label = "") {
  cat(sprintf("\n%s\n", label))
  for (i in seq_len(nrow(res))) {
    cat(sprintf("  %s %s (%s) %s\n",
                if (res$pass[i]) "PASS" else "FAIL",
                res$check[i], res$rule[i],
                if (res$pass[i]) "" else paste("-", res$message[i])))
  }
  invisible(all(res$pass))
}
