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
