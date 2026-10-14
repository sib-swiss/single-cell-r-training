#!/usr/bin/env Rscript
#
# Full validation entrypoint (pack.yml: validity.validation_entrypoint).
# Needs Seurat + edgeR + limma - runs in sibswiss/training-singlecell-rstudio.

here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "de_markers.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

# --- 1. Static lint of the pack's own code -----------------------------------

impl <- file.path(pack, "implementation", "de_markers.R")
lint <- lint_de_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/de_markers.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

# --- 2. Happy path: pseudobulk DE finds the true signal ----------------------

cat("\n=== pseudobulk_de(): two conditions, 3 replicates each ===\n")
seu <- make_de_fixture(c("A", "B"), n_samples_per_cond = 3)
de <- pseudobulk_de(seu, sample_col = "orig.ident", condition_col = "type",
                    contrast = c("B", "A"))

r <- check_pseudobulk_design(de)
if (!report(r, "pseudobulk_de: design post-conditions")) note("pseudobulk design post-conditions failed")

truede_found <- sum(grepl("^TRUEDE", de$gene[de$adj.P.Val < 0.05]))
cat(sprintf("TRUEDE genes called significant: %d / 10\n", truede_found))
if (truede_found < 7) note("pseudobulk_de failed to recover the planted true-DE genes (found %d/10)", truede_found)

noise_fp <- sum(grepl("^NOISE", de$gene[de$adj.P.Val < 0.05]))
cat(sprintf("NOISE genes called significant (false positives): %d / 10\n", noise_fp))
if (noise_fp > 3) note("pseudobulk_de called an implausible number of patient-noise genes significant (%d/10) - check aggregation", noise_fp)

cat(de_methods_sentence(de), "\n")

# --- 3. R-DE-001 regression: reproduce the course's pseudoreplication effect
# in miniature. A naive per-cell test on the SAME data, ignoring the patient
# structure, is expected to call substantially more NOISE genes "significant"
# than pseudobulk did - this is the failure this pack exists to prevent.

cat("\n=== R-DE-001: naive per-cell test must overstate noise-gene significance ===\n")
seu <- Seurat::NormalizeData(seu, verbose = FALSE)
seu <- Seurat::SetIdent(seu, value = "type")
naive <- suppressWarnings(Seurat::FindMarkers(seu, ident.1 = "B", ident.2 = "A",
                                              logfc.threshold = 0, min.pct = 0))
naive_noise_fp <- sum(grepl("^NOISE", rownames(naive)[naive$p_val_adj < 0.05]))
cat(sprintf("NOISE genes called significant by naive per-cell test: %d / 10\n", naive_noise_fp))
cat(sprintf("NOISE genes called significant by pseudobulk:          %d / 10\n", noise_fp))
if (naive_noise_fp <= noise_fp) {
  note("R-DE-001 regression: naive per-cell test did not overstate noise-gene significance relative to pseudobulk on this fixture")
  cat("  FAIL - the fixture/assertion no longer demonstrates the failure mode\n")
} else {
  cat("  PASS - naive per-cell testing overstates significance on patient-level noise, as expected\n")
}

# --- 4. Replicate guard: too few replicates must be refused -----------------

cat("\n=== R-DE-001: too few replicates must be refused ===\n")
seu_sparse <- make_de_fixture(c("A", "B"), n_samples_per_cond = 1)
got <- tryCatch({
  pseudobulk_de(seu_sparse, sample_col = "orig.ident", condition_col = "type",
               contrast = c("B", "A"))
  "no error"
}, error = function(e) conditionMessage(e))
if (identical(got, "no error")) {
  note("pseudobulk_de() accepted only 1 replicate per condition")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required\n")
}

# --- 5. Inconsistent replicate_col mapping must be refused ------------------

cat("\n=== Inconsistent replicate_col mapping must be refused ===\n")
seu_bad <- make_de_fixture(c("A", "B"), n_samples_per_cond = 3)
seu_bad$patient.id <- seu_bad$orig.ident
# Corrupt the mapping: split ONE sample's own cells across two different
# patient IDs, so sample_col no longer maps 1:1 to replicate_col.
first_sample <- seu_bad$orig.ident[1]
idx <- which(seu_bad$orig.ident == first_sample)
half_idx <- idx[seq_len(length(idx) %/% 2)]
seu_bad$patient.id[half_idx] <- "BOGUS-PATIENT"
got2 <- tryCatch({
  pseudobulk_de(seu_bad, sample_col = "orig.ident", condition_col = "type",
               replicate_col = "patient.id", contrast = c("B", "A"))
  "no error"
}, error = function(e) conditionMessage(e))
if (identical(got2, "no error")) {
  note("pseudobulk_de() accepted an inconsistent sample->replicate mapping")
  cat("  FAIL\n")
} else if (!grepl("map to exactly one", got2)) {
  note("pseudobulk_de() errored on the corrupted mapping, but not via the intended structural guard: %s", got2)
  cat("  FAIL (wrong error):", got2, "\n")
} else {
  cat("  PASS - errored as required via the structural mapping guard\n")
}

# --- 6. anova_de() smoke test: 3 groups ------------------------------------

cat("\n=== anova_de(): 3 groups ===\n")
seu3 <- make_de_fixture(c("young", "middle", "old"), n_samples_per_cond = 3)
res3 <- anova_de(seu3, sample_col = "orig.ident", group_col = "type")
if (!("gene" %in% names(res3)) || !nrow(res3)) {
  note("anova_de() did not return a usable result")
} else {
  cat(sprintf("  PASS - anova_de() returned %d genes\n", nrow(res3)))
}

# --- 7. Diagnostics structure (R-DE-005) ------------------------------------

cat("\n=== de_diagnostic_plots() structure ===\n")
plots <- de_diagnostic_plots(de)
rd <- check_de_diagnostics(plots)
if (!report(rd, "diagnostic plot post-conditions")) note("de_diagnostic_plots() missing expected components")

# --- verdict -----------------------------------------------------------------

cat("\n", strrep("-", 60), "\n", sep = "")
if (length(failures)) {
  cat("VALIDATION FAILED\n")
  for (f in failures) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("VALIDATION PASSED\n")
cat("Validated against: R ", as.character(getRversion()),
    ", Seurat ", as.character(packageVersion("Seurat")),
    ", edgeR ", as.character(packageVersion("edgeR")),
    ", limma ", as.character(packageVersion("limma")), "\n", sep = "")
quit(status = 0)
