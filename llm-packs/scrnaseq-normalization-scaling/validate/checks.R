#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on a scale_variable_features() result
check_scaled_object <- function(seu) {
  res <- list()

  has_data <- !is.null(SeuratObject::LayerData(seu, layer = "data", assay = SeuratObject::DefaultAssay(seu)))
  res[[length(res)+1]] <- chk("C-NORM-001", "manual", has_data, "no normalised 'data' layer found")

  vf <- Seurat::VariableFeatures(seu)
  res[[length(res)+1]] <- chk("C-NORM-002", "manual", length(vf) > 0,
                              "no variable features recorded")

  scale_dat <- tryCatch(SeuratObject::LayerData(seu, layer = "scale.data",
                                                assay = SeuratObject::DefaultAssay(seu)),
                        error = function(e) NULL)
  res[[length(res)+1]] <- chk("C-NORM-003", "manual",
                              !is.null(scale_dat) && all(vf %in% rownames(scale_dat)),
                              "scale.data does not cover all variable features")

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
