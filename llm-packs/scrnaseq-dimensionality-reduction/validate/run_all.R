#!/usr/bin/env Rscript
here <- dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
pack <- normalizePath(file.path(here, ".."))

source(file.path(here, "checks.R"))
source(file.path(pack, "implementation", "dimred.R"))
source(file.path(here, "make_fixtures.R"))
source(file.path(here, "lint_code.R"))

failures <- character()
note <- function(...) failures <<- c(failures, sprintf(...))

impl <- file.path(pack, "implementation", "dimred.R")
lint <- lint_dimred_code(readLines(impl, warn = FALSE))
cat(format_findings(lint, "implementation/dimred.R"))
if (any(lint$severity == "block")) note("reference implementation has blocking lint findings")

cat("\n=== Happy path: run_dimred() ===\n")
seu <- make_fixture()
seu_d <- run_dimred(seu, dims = 1:20)
r <- check_reduction_object(seu_d)
if (!report(r, "run_dimred(): post-conditions")) note("post-conditions failed")
cat(dimred_methods_sentence(seu_d), "\n")

cat("\n=== R-DR-002: refuses to run without a PCA reduction ===\n")
seu_nopca <- make_fixture(run_pca = FALSE)
got <- tryCatch({ run_dimred(seu_nopca, dims = 1:20); "no error" },
                error = function(e) conditionMessage(e))
if (identical(got, "no error")) {
  note("R-DR-002 REGRESSION: run_dimred() proceeded with no PCA reduction present")
  cat("  FAIL - this is the failure the pack exists to prevent\n")
} else {
  cat("  PASS - errored as required:", got, "\n")
}

cat("\n=== dims required ===\n")
got2 <- tryCatch({ run_dimred(seu); "no error" }, error = function(e) conditionMessage(e))
if (identical(got2, "no error")) {
  note("run_dimred() accepted a missing dims argument")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required\n")
}

cat("\n=== dims cannot exceed computed PCs ===\n")
got3 <- tryCatch({ run_dimred(seu, dims = 1:100); "no error" },
                 error = function(e) conditionMessage(e))
if (identical(got3, "no error")) {
  note("run_dimred() accepted dims beyond the number of computed PCs")
  cat("  FAIL\n")
} else {
  cat("  PASS - errored as required:", got3, "\n")
}

cat("\n=== dims honoured on a differently-named (integrated) reduction ===\n")
seu_int <- seu
seu_int@reductions[["integrated.cca"]] <- seu_int@reductions[["pca"]]
seu_i <- run_dimred(seu_int, dims = 1:20, reduction = "integrated.cca")
if (!"umap" %in% names(seu_i@reductions)) {
  note("run_dimred() did not produce a UMAP when using a named integrated reduction")
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
