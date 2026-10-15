#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on a cluster_cells() result
check_clustered_object <- function(seu) {
  res <- list()

  prov <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  integrated <- grep("integrated|harmony|cca|mnn", names(seu@reductions),
                     ignore.case = TRUE, value = TRUE)
  res[[length(res)+1]] <- chk(
    "C-CLU-002", "R-CLU-002",
    !length(integrated) || (!is.null(prov) && prov$reduction %in% integrated),
    "object carries an integrated reduction but clustering did not use it")

  cluster_cols <- grep("_snn_res\\.", colnames(seu@meta.data))
  res[[length(res)+1]] <- chk("C-CLU-010", "manual", length(cluster_cols) > 0,
                              "no *_snn_res.* columns found - did FindClusters() run?")

  do.call(rbind, res)
}

#' Structural post-conditions on a check_sample_mixing() result
check_sample_mixing_result <- function(res) {
  r <- list()

  r[[length(r)+1]] <- chk(
    "C-CLU-006a", "R-CLU-006",
    is.data.frame(res$purity) &&
      all(c("cluster", "n_cells", "purity", "dominant_sample") %in% names(res$purity)),
    "purity data.frame missing expected columns")

  r[[length(r)+1]] <- chk(
    "C-CLU-006b", "R-CLU-006",
    is.data.frame(res$purity) && all(res$purity$purity >= 0 & res$purity$purity <= 1),
    "purity values out of the [0, 1] range")

  r[[length(r)+1]] <- chk(
    "C-CLU-006c", "R-CLU-006", inherits(res$centroid_distances, "dist"),
    "centroid_distances is not a dist object")

  do.call(rbind, r)
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
