#!/usr/bin/env Rscript
#
# Full validation entrypoint (pack.yml: validity.validation_entrypoint).
# Needs Seurat - runs in sibswiss/training-singlecell-rstudio.
#
# On success, CI stamps validity.last_validated in pack.yml. That stamp is the
# pack's whole claim to trustworthiness, so nothing else may write it.

here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "qc_filter.R"))
source(file.path(here, "make_fixtures.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

# --- 1. Static lint of the pack's own code -----------------------------------
# The reference implementation must satisfy the rules it teaches.

source(file.path(here, "lint_code.R"))
impl <- file.path(pack, "implementation", "qc_filter.R")
lint <- lint_qc_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/qc_filter.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

# --- 2. Happy path, three species --------------------------------------------

for (sp in c("human", "mouse", "drosophila")) {
  cat(sprintf("\n=== %s ===\n", sp))
  seu <- make_fixture(sp)
  seu <- add_qc_metrics(seu, species = sp)

  r <- check_qc_object(seu)
  if (!report(r, sprintf("%s: QC metric post-conditions", sp))) {
    note("%s: QC post-conditions failed", sp)
  }

  print(qc_summary(seu))

  thr <- list(min_features = 50, max_features = 5000, max_mito = 25)
  seu_f <- filter_cells(seu, thr$min_features, thr$max_features, thr$max_mito)
  rf <- check_filtered_object(seu, seu_f)
  if (!report(rf, sprintf("%s: filtering post-conditions", sp))) {
    note("%s: filtering post-conditions failed", sp)
  }
}

# --- 3. The silent-failure mode must be LOUD ---------------------------------
# This is the test the whole pack exists for: applying the human pattern to
# mouse data must stop, not sail through with percent.mito == 0.

cat("\n=== R-001: wrong species must error, not pass silently ===\n")
mouse <- make_fixture("mouse")
got <- tryCatch({ add_qc_metrics(mouse, species = "human"); "no error" },
                error = function(e) conditionMessage(e))
if (identical(got, "no error")) {
  note("R-001 REGRESSION: human pattern on mouse data did not error")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required\n")
}

cat("\n=== R-001: species argument must be required ===\n")
got2 <- tryCatch({ add_qc_metrics(mouse); "no error" },
                 error = function(e) conditionMessage(e))
if (identical(got2, "no error")) {
  note("add_qc_metrics() accepted a missing species argument")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required\n")
}

# --- 4. R-013: filtering after integration must be refused -------------------

cat("\n=== R-013: filtering post-integration must be refused ===\n")
h <- add_qc_metrics(make_fixture("human"), species = "human")
h@reductions$integrated.cca <- SeuratObject::CreateDimReducObject(
  embeddings = matrix(0, ncol(h), 2,
                      dimnames = list(colnames(h), c("cca_1","cca_2"))),
  key = "cca_", assay = "RNA")
got3 <- tryCatch({ filter_cells(h, 50, 5000, 25); "no error" },
                 error = function(e) conditionMessage(e))
if (identical(got3, "no error")) {
  note("R-013: filtering an integrated object was allowed")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required\n")
}

# --- 5. R-009: QC metrics computed after normalisation must be refused ------

cat("\n=== R-009: add_qc_metrics() on an already-normalised object must error ===\n")
norm <- make_fixture("human")
norm <- Seurat::NormalizeData(norm, verbose = FALSE)
got4 <- tryCatch({ add_qc_metrics(norm, species = "human"); "no error" },
                 error = function(e) conditionMessage(e))
if (identical(got4, "no error")) {
  note("R-009 REGRESSION: add_qc_metrics() did not notice the object was already normalised")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required:", got4, "\n")
}

# --- 6. R-004/P-012: filter_cells() on unjoined multi-sample layers ---------

cat("\n=== R-004/P-012: filter_cells() on unjoined layers must error, not silently cover one sample ===\n")
unjoined <- add_qc_metrics(make_fixture("human"), species = "human")
unjoined[["RNA"]] <- split(unjoined[["RNA"]], f = unjoined$orig.ident)
got5 <- tryCatch({ filter_cells(unjoined, 50, 5000, 25); "no error" },
                 error = function(e) conditionMessage(e))
if (identical(got5, "no error")) {
  note("R-004/P-012 REGRESSION: filter_cells() accepted an object with unjoined per-sample layers")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required:", got5, "\n")
}

# --- verdict -----------------------------------------------------------------

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
cat("Validated against: R ", as.character(getRversion()),
    ", Seurat ", as.character(packageVersion("Seurat")), "\n", sep = "")
quit(status = 0)
