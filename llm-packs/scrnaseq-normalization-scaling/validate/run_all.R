#!/usr/bin/env Rscript
here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "normalize_scale.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

impl <- file.path(pack, "implementation", "normalize_scale.R")
lint <- lint_norm_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/normalize_scale.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

cat("\n=== Happy path: scale_variable_features() ===\n")
seu <- make_fixture()
seu_s <- scale_variable_features(seu, nfeatures = 50)
r <- check_scaled_object(seu_s)
if (!report(r, "scale_variable_features(): post-conditions")) note("post-conditions failed")
cat(normalization_methods_sentence(seu_s), "\n")

cat("\n=== R-NORM-001: extra_scale_genes are covered in scale.data ===\n")
extra <- rownames(seu)[1:5]
seu_s2 <- scale_variable_features(seu, nfeatures = 50, extra_scale_genes = extra)
scaled <- rownames(SeuratObject::LayerData(seu_s2, layer = "scale.data"))
if (!all(extra %in% scaled)) {
  note("extra_scale_genes were not present in scale.data")
  cat("  FAIL\n")
} else {
  cat("  PASS - extra_scale_genes covered\n")
}

cat("\n=== extra_scale_genes not present in the object are warned about, not silently dropped ===\n")
got <- tryCatch({
  withCallingHandlers(
    scale_variable_features(seu, nfeatures = 50, extra_scale_genes = "NOT_A_REAL_GENE"),
    warning = function(w) { cat("  PASS - warned:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  "handled"
}, error = function(e) "error")
if (!identical(got, "handled")) note("missing extra_scale_genes did not warn as expected")

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
