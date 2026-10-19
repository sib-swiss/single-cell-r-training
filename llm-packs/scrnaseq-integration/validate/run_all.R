#!/usr/bin/env Rscript
#
# Full validation entrypoint (pack.yml: validity.validation_entrypoint).
# Needs Seurat - runs in sibswiss/training-singlecell-rstudio.
#
# On success, CI stamps validity.last_validated in pack.yml.

here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "integrate.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

# --- 1. Static lint of the pack's own code -----------------------------------

impl <- file.path(pack, "implementation", "integrate.R")
lint <- lint_integration_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/integrate.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

# --- 2. Happy path: RNA-assay integration, condition crossed with batch ------

cat("\n=== RNA-assay integration ===\n")
seu <- make_integration_fixture()
seu <- Seurat::NormalizeData(seu, verbose = FALSE)
seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
seu <- Seurat::ScaleData(seu, verbose = FALSE)
seu <- Seurat::RunPCA(seu, npcs = 10, verbose = FALSE)

seu_i <- integrate_samples(seu, assay = "RNA", batch_var = "orig.ident")

r <- check_integrated_object(seu_i, assay = "RNA", reduction = "integrated.cca")
if (!report(r, "RNA integration: post-conditions")) note("RNA integration post-conditions failed")

cat(integration_methods_sentence(seu_i), "\n")

# --- 3. R-INT-003: the planted biological signal must still separate --------
# condition was crossed with batch (not confounded) when the fixture was
# built, specifically so this positive control has something real to find.

cat("\n=== R-INT-003: condition signal must survive integration ===\n")
sig <- check_integration_preserves_signal(seu_i, reduction = "integrated.cca",
                                          group_by = "condition")
print(sig)
if (!any(sig$separates, na.rm = TRUE)) {
  note("planted condition signal did not survive integration - over-integration regression")
  cat("  FAIL - this is the failure the diagnostic exists to catch\n")
} else {
  cat("  PASS - condition still separates post-integration\n")
}

# --- 4. R-INT-005: wrong-assay PCA must be refused, not silently used -------

cat("\n=== R-INT-005: SCT integration against an RNA-assay PCA must error ===\n")
seu2 <- make_integration_fixture()
seu2 <- Seurat::NormalizeData(seu2, verbose = FALSE)
seu2 <- Seurat::FindVariableFeatures(seu2, verbose = FALSE)
seu2 <- Seurat::ScaleData(seu2, verbose = FALSE)
seu2 <- Seurat::RunPCA(seu2, npcs = 10, verbose = FALSE)   # "pca", assay.used = RNA
seu2[["SCT"]] <- seu2[["RNA"]]                             # stand-in SCT assay for the test

got <- tryCatch({
  integrate_samples(seu2, assay = "SCT", batch_var = "orig.ident", orig_reduction = "pca")
  "no error"
}, error = function(e) conditionMessage(e))
if (identical(got, "no error")) {
  note("R-INT-005 REGRESSION: SCT integration against an RNA-assay PCA did not error")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required:", got, "\n")
}

# --- 4b. SCT-assay happy path against a REAL SCTransform() output -----------
# The test above only exercises the wrong-assay guard with a stand-in SCT
# assay (a copy of RNA). A genuine SCTransform() output is a different Seurat
# class (SCTAssay) that split() cannot handle directly without first
# converting it to a plain Assay5 - a real bug found by hand-testing this
# pack's own SKILL.md example, not by this fixture, until now.

cat("\n=== SCT-assay integration against a REAL SCTransform() output ===\n")
seu3 <- make_integration_fixture()
seu3 <- Seurat::NormalizeData(seu3, verbose = FALSE)
seu3 <- Seurat::FindVariableFeatures(seu3, verbose = FALSE)
seu3 <- Seurat::ScaleData(seu3, verbose = FALSE)
seu3 <- Seurat::RunPCA(seu3, npcs = 10, verbose = FALSE)
seu3 <- Seurat::SCTransform(seu3, verbose = FALSE)
Seurat::DefaultAssay(seu3) <- "SCT"
seu3 <- Seurat::RunPCA(seu3, reduction.name = "pcaSCT", verbose = FALSE)

got2 <- tryCatch({
  seu3_i <- integrate_samples(seu3, assay = "SCT", batch_var = "orig.ident",
                              orig_reduction = "pcaSCT", new_reduction = "SCT.cca")
  if (!"SCT.cca" %in% names(seu3_i@reductions)) "no SCT.cca reduction produced" else "ok"
}, error = function(e) conditionMessage(e))
if (!identical(got2, "ok")) {
  note("SCT-assay integration against a real SCTransform() output failed: %s", got2)
  cat("  FAIL:", got2, "\n")
} else {
  cat("  PASS\n")
}

# --- 5. Single-batch object must be refused ---------------------------------

cat("\n=== Single-batch object must be refused ===\n")
seu3 <- make_integration_fixture()
seu3 <- seu3[, seu3$orig.ident == "S1"]
seu3 <- Seurat::NormalizeData(seu3, verbose = FALSE)
seu3 <- Seurat::FindVariableFeatures(seu3, verbose = FALSE)
seu3 <- Seurat::ScaleData(seu3, verbose = FALSE)
seu3 <- Seurat::RunPCA(seu3, npcs = 10, verbose = FALSE)
got2 <- tryCatch({
  integrate_samples(seu3, assay = "RNA", batch_var = "orig.ident")
  "no error"
}, error = function(e) conditionMessage(e))
if (identical(got2, "no error")) {
  note("integrate_samples() accepted a single-batch object")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required\n")
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
