#!/usr/bin/env Rscript
#
# Runtime post-conditions from rules/integration.rules.yml. Needs Seurat, so
# this runs in the course container, not on a bare R image.

suppressPackageStartupMessages(library(Seurat))

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on an object that has had integrate_samples() applied
check_integrated_object <- function(seu, assay, reduction) {
  res <- list()

  # R-INT-004: layers must be joined, not left per-sample.
  layers <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) character(0))
  n_count_layers <- sum(grepl("^counts", layers))
  res[[length(res)+1]] <- chk(
    "C-INT-004", "R-INT-004", n_count_layers <= 1,
    sprintf("%d counts layers present on assay '%s'; JoinLayers() should leave exactly one",
            n_count_layers, assay))

  # The integrated reduction must actually exist and have >0 dimensions.
  has_reduction <- reduction %in% names(seu@reductions)
  res[[length(res)+1]] <- chk(
    "C-INT-010", "manual", has_reduction,
    sprintf("expected reduction '%s' not found on the object", reduction))

  if (has_reduction) {
    ndims <- ncol(SeuratObject::Embeddings(seu, reduction))
    res[[length(res)+1]] <- chk(
      "C-INT-011", "manual", ndims > 0,
      sprintf("reduction '%s' has no dimensions", reduction))
  }

  # R-INT-001: the default assay must not be left on a legacy "integrated"
  # assay produced by IntegrateData() - that assay is for clustering/
  # visualisation only and must not look like the analysis's resting state.
  res[[length(res)+1]] <- chk(
    "C-INT-001", "R-INT-001",
    !identical(tolower(Seurat::DefaultAssay(seu)), "integrated"),
    "DefaultAssay is 'integrated' after integration; reset to the original assay before returning the object, so DE is not accidentally run on corrected data")

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
