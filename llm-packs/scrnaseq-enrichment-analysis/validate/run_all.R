#!/usr/bin/env Rscript
here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "enrichment.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

impl <- file.path(pack, "implementation", "enrichment.R")
lint <- lint_enrichment_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/enrichment.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

cat("\n=== Happy path: run_ora() finds the planted cell-cycle signal ===\n")
fx <- make_fixture()
ora <- run_ora(fx$genes, universe = fx$universe, ont = "BP", min_gs_size = 5)
r <- check_enrichment_result(ora)
if (!report(r, "run_ora(): post-conditions")) note("post-conditions failed")

sig_terms <- ora@result$Description[ora@result$p.adjust < 0.05]
cat(sprintf("Significant GO terms found: %d\n", length(sig_terms)))
cell_cycle_hit <- any(grepl("cell cycle|mitotic|cycle", sig_terms, ignore.case = TRUE))
if (!cell_cycle_hit) {
  note("run_ora() did not recover any cell-cycle-related GO term for a planted cell-cycle gene set")
  cat("  FAIL\n")
} else {
  cat("  PASS - recovered a cell-cycle-related term\n")
}
cat(enrichment_methods_sentence(ora), "\n")

cat("\n=== R-ENR-001: universe is required, not optional ===\n")
got <- tryCatch({ run_ora(fx$genes); "no error" }, error = function(e) conditionMessage(e))
if (identical(got, "no error")) {
  note("R-ENR-001 REGRESSION: run_ora() accepted a missing universe argument")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required\n")
}

cat("\n=== R-ENR-003: GSEA refuses an unsorted ranked list ===\n")
ranked <- make_ranked_fixture()
shuffled <- sample(ranked)
got2 <- tryCatch({ run_gsea(shuffled); "no error" }, error = function(e) conditionMessage(e))
if (identical(got2, "no error")) {
  note("R-ENR-003 REGRESSION: run_gsea() accepted an unsorted ranked list")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required:", got2, "\n")
}

cat("\n=== run_gsea() on a properly sorted, full ranked list ===\n")
gsea_res <- run_gsea(ranked, ont = "BP", min_full_list_size = 500)
if (is.null(gsea_res) || !is.data.frame(gsea_res@result)) {
  note("run_gsea() did not return a usable result")
} else {
  cat(sprintf("  PASS - %d terms in GSEA result\n", nrow(gsea_res@result)))
}

cat("\n=== run_gsea() warns on a suspiciously short ranked list ===\n")
short_ranked <- ranked[1:50]
got3 <- tryCatch({
  withCallingHandlers(
    run_gsea(short_ranked, min_full_list_size = 1000),
    warning = function(w) { cat("  PASS - warned:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  "handled"
}, error = function(e) "error")
if (!identical(got3, "handled")) note("run_gsea() did not warn on a suspiciously short ranked list")

cat("\n=== append_methods(): creates, then replaces (not duplicates) on re-run ===\n")
methods_path <- tempfile(fileext = ".md")
append_methods("First sentence.", path = methods_path, section = "test:section")
append_methods("Updated sentence.", path = methods_path, section = "test:section")
methods_lines <- readLines(methods_path)
n_blocks <- sum(grepl("^<!-- methods:test:section:begin -->$", methods_lines))
if (n_blocks != 1 || !any(grepl("Updated sentence\\.", methods_lines)) ||
    any(grepl("First sentence\\.", methods_lines))) {
  note("append_methods() duplicated or failed to replace an existing section on re-run")
  cat("  FAIL\n")
} else {
  cat("  PASS - section replaced in place, not duplicated\n")
}
unlink(methods_path)

cat("\n", strrep("-", 60), "\n", sep = "")
if (length(failures)) {
  cat("VALIDATION FAILED\n")
  for (f in failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("VALIDATION PASSED\n")
cat("Validated against: R ", as.character(getRversion()), "\n", sep = "")
quit(status = 0)
