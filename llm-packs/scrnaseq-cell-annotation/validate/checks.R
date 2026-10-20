#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on an annotate_with_reference() result
check_annotation_object <- function(seu) {
  res <- list()

  res[[length(res)+1]] <- chk("C-ANN-001", "manual",
                              !is.null(seu$SingleR_annot),
                              "no SingleR_annot column found in meta.data")

  res[[length(res)+1]] <- chk("C-ANN-002", "R-ANN-002",
                              !is.null(attr(seu, "annotation_result")),
                              "no annotation_result attribute recorded - diagnostics cannot be generated")

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
