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

cat("\n=== R-ANN-007: match_clusters_to_markers() recovers the planted cell type ===\n")
# Synthetic marker_db shaped like a real CellMarker/PanglaoDB table (a
# cell-type column + a gene column), built from the fixture's own planted
# ground truth plus one decoy cell type with unrelated genes - this stays
# testable in CI without a network download.
cluster_markers <- data.frame(
  cluster = rep(names(fx$marker_idx), each = 20),
  gene = fx$genes[unlist(fx$marker_idx)],
  stringsAsFactors = FALSE)
decoy_genes <- setdiff(fx$genes, fx$genes[unlist(fx$marker_idx)])[1:10]
marker_db <- data.frame(
  cell_name = c(rep(names(fx$marker_idx), each = 20), rep("Unrelated_decoy", 10)),
  Symbol = c(fx$genes[unlist(fx$marker_idx)], decoy_genes),
  stringsAsFactors = FALSE)

matches <- match_clusters_to_markers(cluster_markers, marker_db,
                                     celltype_col = "cell_name", gene_col = "Symbol")
print(matches)
if (!all(matches$best_match == matches$cluster)) {
  note("match_clusters_to_markers() failed to recover the planted cell type for at least one cluster")
  cat("  FAIL\n")
} else {
  cat("  PASS - recovered the planted cell type for every cluster\n")
}

cat("\n=== match_clusters_to_markers(): errors on a missing column name ===\n")
got_mm <- tryCatch({
  match_clusters_to_markers(cluster_markers, marker_db,
                            celltype_col = "nonexistent_col", gene_col = "Symbol")
  "no error"
}, error = function(e) conditionMessage(e))
if (identical(got_mm, "no error")) {
  note("match_clusters_to_markers() accepted a nonexistent celltype_col")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required:", got_mm, "\n")
}

cat("\n=== stamp_provenance()/annotation_methods_sentence(): marker_source path ===\n")
seu_manual <- stamp_provenance(fx$seu, marker_source = "CellMarker")
cat(annotation_methods_sentence(seu_manual), "\n")
seu_memory <- stamp_provenance(fx$seu, marker_source = "model_memory_unverified")
sentence_memory <- annotation_methods_sentence(seu_memory)
cat(sentence_memory, "\n")
if (!grepl("NOT independently verified", sentence_memory)) {
  note("annotation_methods_sentence() did not flag the model_memory_unverified fallback as unverified")
}

cat("\n=== Rare labels get relabelled to 'none', not silently dropped ===\n")
# Force a near-empty label by using a very high min_label_cells.
seu_rare <- annotate_with_reference(fx$seu, ref = fx$ref, labels = fx$labels, min_label_cells = 1000)
if (!all(seu_rare$SingleR_annot == "none")) {
  note("min_label_cells did not relabel everything when set above the total cell count")
} else {
  cat("  PASS\n")
}

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
