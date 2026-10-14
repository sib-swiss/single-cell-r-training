#!/usr/bin/env Rscript
chk <- function(id, rule, ok, msg) {
  data.frame(check = id, rule = rule, pass = isTRUE(ok), message = msg,
             stringsAsFactors = FALSE)
}

#' Post-conditions on a run_ora() result
check_enrichment_result <- function(res) {
  rs <- list()
  prov <- utils::tail(attr(res, "provenance"), 1)[[1]]

  rs[[length(rs)+1]] <- chk("C-ENR-001", "R-ENR-001",
                            !is.null(prov) && !is.null(prov$n_universe) && prov$n_universe > 0,
                            "no universe recorded in provenance - was this run through run_ora()?")

  rs[[length(rs)+1]] <- chk("C-ENR-010", "manual",
                            !is.null(res) && is.data.frame(res@result),
                            "no @result data frame found on the enrichment object")

  do.call(rbind, rs)
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
