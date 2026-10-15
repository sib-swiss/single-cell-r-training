#!/usr/bin/env Rscript
here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "cluster.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

impl <- file.path(pack, "implementation", "cluster.R")
lint <- lint_clustering_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/cluster.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

cat("\n=== Happy path: single-sample object (no integration) ===\n")
seu <- make_fixture(with_integration = FALSE)
seu_c <- cluster_cells(seu, dims = 1:15)
r <- check_clustered_object(seu_c)
if (!report(r, "cluster_cells(): post-conditions")) note("post-conditions failed")
cat(clustering_methods_sentence(seu_c), "\n")

cat("\n=== R-CLU-002: refuses to default to pca on an integrated object ===\n")
seu_int <- make_fixture(with_integration = TRUE)
got <- tryCatch({ cluster_cells(seu_int, dims = 1:15); "no error" },
                error = function(e) conditionMessage(e))
if (identical(got, "no error")) {
  note("R-CLU-002 REGRESSION: cluster_cells() silently defaulted to pca on an integrated object")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required:", got, "\n")
}

cat("\n=== Explicit integrated reduction is honoured and recorded ===\n")
seu_int_c <- cluster_cells(seu_int, dims = 1:15, reduction = "integrated.cca")
r2 <- check_clustered_object(seu_int_c)
if (!report(r2, "cluster_cells() on integrated reduction: post-conditions")) {
  note("integrated-reduction post-conditions failed")
}

cat("\n=== clustree helper runs on a resolution sweep ===\n")
got2 <- tryCatch({ p <- cluster_stability_tree(seu_c); "ok" },
                 error = function(e) conditionMessage(e))
if (!identical(got2, "ok")) {
  note("cluster_stability_tree() failed: %s", got2)
  cat("  FAIL:", got2, "\n")
} else {
  cat("  PASS\n")
}

# --- R-CLU-006: check_sample_mixing() must discriminate a planted sample
# effect from a no-effect fixture, on the actual UMAP embedding -----------

cat("\n=== R-CLU-006: check_sample_mixing() flags a planted sample/batch effect ===\n")
seu_mix <- make_fixture(with_integration = TRUE)
seu_mix <- cluster_cells(seu_mix, dims = 1:15, reduction = "integrated.cca")
mix <- check_sample_mixing(seu_mix, sample_col = "orig.ident")
r3 <- check_sample_mixing_result(mix)
if (!report(r3, "check_sample_mixing(): structural post-conditions")) {
  note("check_sample_mixing() structure check failed")
}
mean_purity_effect <- mean(mix$purity$purity)
cat(sprintf("  mean purity (planted sample effect): %.2f\n", mean_purity_effect))
if (mean_purity_effect < 0.7) {
  note("check_sample_mixing() did not detect the planted sample/batch effect via purity (mean purity %.2f)", mean_purity_effect)
}

cat("\n=== check_sample_mixing(): stays low with no planted sample effect ===\n")
seu_nomix <- make_fixture(with_integration = FALSE)
seu_nomix <- cluster_cells(seu_nomix, dims = 1:15)
mix2 <- check_sample_mixing(seu_nomix, sample_col = "orig.ident")
mean_purity_nomix <- mean(mix2$purity$purity)
cat(sprintf("  mean purity (no sample effect):      %.2f\n", mean_purity_nomix))
if (mean_purity_nomix >= mean_purity_effect) {
  note("check_sample_mixing(): no-effect fixture did not have lower purity than the planted-effect fixture (%.2f vs %.2f) - diagnostic not discriminating", mean_purity_nomix, mean_purity_effect)
} else {
  cat("  PASS - planted-effect fixture shows higher purity than the no-effect fixture\n")
}

cat("\n=== check_sample_mixing(): refuses a single-sample object ===\n")
seu_single <- make_fixture(with_integration = FALSE)
seu_single <- seu_single[, seu_single$orig.ident == "S1"]
seu_single <- cluster_cells(seu_single, dims = 1:15)
got4 <- tryCatch({ check_sample_mixing(seu_single); "no error" },
                 error = function(e) conditionMessage(e))
if (identical(got4, "no error")) {
  note("check_sample_mixing() accepted a single-sample object")
  cat("  FAIL\n")
} else if (!grepl("Only one sample/batch", got4)) {
  note("check_sample_mixing() errored on a single-sample object, but not via the intended guard: %s", got4)
  cat("  FAIL (wrong error):", got4, "\n")
} else {
  cat("  PASS - errored as required via the sample-count guard:", got4, "\n")
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
