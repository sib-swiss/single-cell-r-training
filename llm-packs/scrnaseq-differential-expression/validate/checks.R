#!/usr/bin/env Rscript
#
# Runtime post-conditions from rules/differential-expression.rules.yml.
# Operates on the data.frame results returned by implementation/de_markers.R,
# so (unlike the QC and Integration packs' checks.R) this one does not itself
# require Seurat to be loaded - only whatever produced the result did.

chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on a pseudobulk_de() result
check_pseudobulk_design <- function(de_result, min_replicates = 2) {
  res <- list()
  prov <- utils::tail(attr(de_result, "provenance"), 1)[[1]]

  res[[length(res)+1]] <- chk(
    "C-DE-001", "manual", !is.null(prov) && identical(prov$step, "pseudobulk_de"),
    "no pseudobulk_de provenance found on this result - was it produced by pseudobulk_de()?")

  if (!is.null(prov) && !is.null(prov$n_replicates)) {
    enough <- all(unlist(prov$n_replicates) >= min_replicates)
    res[[length(res)+1]] <- chk(
      "C-DE-002", "R-DE-001", enough,
      sprintf("fewer than %d replicates in at least one condition (R-DE-001 pseudoreplication guard)",
              min_replicates))
  }

  # R-DE-006 is a reporting-style rule (manual only - see rules.yml) and is
  # not checked here; nothing here can verify ident.1/ident.2 were stated in
  # the user-facing reply.

  do.call(rbind, res)
}

#' Structural check that diagnostics were actually generated (R-DE-005)
check_de_diagnostics <- function(plots) {
  res <- list()
  res[[length(res)+1]] <- chk(
    "C-DE-010", "R-DE-005",
    is.list(plots) && !is.null(plots$p_value_histogram),
    "p-value histogram missing from de_diagnostic_plots() output")
  res[[length(res)+1]] <- chk(
    "C-DE-011", "R-DE-005",
    is.list(plots) && "ma_plot" %in% names(plots),
    "MA plot slot missing from de_diagnostic_plots() output")
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
