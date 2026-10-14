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
