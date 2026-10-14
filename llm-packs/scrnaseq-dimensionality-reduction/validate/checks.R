#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on a run_dimred() result
check_reduction_object <- function(seu) {
  res <- list()

  res[[length(res)+1]] <- chk("C-DR-002", "R-DR-002", "pca" %in% names(seu@reductions),
                              "expected a 'pca' reduction to still be present")
  res[[length(res)+1]] <- chk("C-DR-010", "manual", "umap" %in% names(seu@reductions),
                              "expected a 'umap' reduction after run_dimred()")

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
