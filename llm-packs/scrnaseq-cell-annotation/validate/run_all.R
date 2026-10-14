#!/usr/bin/env Rscript
here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "annotate.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

impl <- file.path(pack, "implementation", "annotate.R")
lint <- lint_annotation_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/annotate.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

cat("\n=== Happy path: annotate_with_reference() ===\n")
fx <- make_fixture()
seu_a <- annotate_with_reference(fx$seu, ref = fx$ref, labels = fx$labels, min_label_cells = 5)
r <- check_annotation_object(seu_a)
if (!report(r, "annotate_with_reference(): post-conditions")) note("post-conditions failed")

accuracy <- mean(seu_a$SingleR_annot == seu_a$true_type)
cat(sprintf("Label accuracy against planted true_type: %.1f%%\n", 100 * accuracy))
if (accuracy < 0.8) note("annotation accuracy on an easy synthetic fixture was implausibly low (%.1f%%)", 100 * accuracy)

cat(annotation_methods_sentence(seu_a), "\n")

cat("\n=== R-ANN-002: no rejection option - decoy reference type is never used when absent from the query ===\n")
if ("NK_cell_DECOY" %in% seu_a$SingleR_annot) {
  note("decoy reference type was assigned despite no decoy cells in the query - unexpected for this fixture")
} else {
  cat("  PASS - no cell matched the decoy type on this easy fixture, as expected\n")
}
cat("  (Every cell still receives a best-available label even when its true type is absent from a\n")
cat("   reference - SingleR itself documents no rejection option; this fixture cannot force that\n")
cat("   case without a genuinely out-of-reference cell population, which is exactly the hazard R-ANN-002 names.)\n")

cat("\n=== Diagnostic plots run without erroring ===\n")
got <- tryCatch({ p <- annotation_diagnostic_plots(seu_a); "ok" },
                error = function(e) conditionMessage(e))
if (!identical(got, "ok")) {
  note("annotation_diagnostic_plots() failed: %s", got)
  cat("  FAIL:", got, "\n")
} else {
  cat("  PASS\n")
}

cat("\n=== Rare labels get relabelled to 'none', not silently dropped ===\n")
# Force a near-empty label by using a very high min_label_cells.
seu_rare <- annotate_with_reference(fx$seu, ref = fx$ref, labels = fx$labels, min_label_cells = 1000)
if (!all(seu_rare$SingleR_annot == "none")) {
  note("min_label_cells did not relabel everything when set above the total cell count")
} else {
  cat("  PASS\n")
}

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
