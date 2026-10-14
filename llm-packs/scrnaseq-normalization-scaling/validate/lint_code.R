#!/usr/bin/env Rscript
#
# Static checks for scrnaseq-normalization-scaling.

strip_comments <- function(lines) {
  vapply(lines, function(l) {
    chars <- strsplit(l, "", fixed = TRUE)[[1]]
    if (!length(chars)) return("")
    in_str <- FALSE; quote_ch <- ""
    for (i in seq_along(chars)) {
      ch <- chars[i]
      if (in_str) {
        if (ch == quote_ch && (i == 1 || chars[i - 1] != "\\")) in_str <- FALSE
      } else if (ch %in% c("'", '"')) {
        in_str <- TRUE; quote_ch <- ch
      } else if (ch == "#") {
        return(paste(chars[seq_len(i - 1)], collapse = ""))
      }
    }
    l
  }, character(1), USE.NAMES = FALSE)
}

finding <- function(rule, lint, severity, line, message) {
  data.frame(rule = rule, lint = lint, severity = severity,
             line = line, message = message, stringsAsFactors = FALSE)
}

no_findings <- function() {
  data.frame(rule = character(), lint = character(), severity = character(),
             line = integer(), message = character(), stringsAsFactors = FALSE)
}

lint_norm_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  blob  <- paste(code, collapse = "\n")
  out   <- list()

  # L-NORM-001 (R-NORM-001): bare ScaleData(seu) with no features=, in a file
  # that also calls DoHeatmap() - a concrete signal that scaled values for a
  # specific gene panel will be needed, which a variable-features-only scale
  # will not have unless that panel happens to be fully inside it.
  scale_lines <- grep("ScaleData\\s*\\(", code)
  has_heatmap <- grepl("DoHeatmap\\s*\\(", blob)
  for (i in scale_lines) {
    bare <- !grepl("features\\s*=", code[i])
    if (bare && has_heatmap) {
      out[[length(out) + 1]] <- finding(
        "R-NORM-001", "L-NORM-001", "warn", i,
        "ScaleData() with no features= argument only scales variable features, but this file also calls DoHeatmap() - if the heatmap's genes aren't all variable features, pass features= explicitly (e.g. union(VariableFeatures(seu), my_genes)).")
    }
  }

  # L-NORM-002 (R-NORM-004): SCTransform() with no later DefaultAssay reset.
  sct_lines <- grep("SCTransform\\s*\\(", code)
  reset_lines <- grep('DefaultAssay\\s*\\([^)]*\\)\\s*<-\\s*["\']RNA["\']', code)
  for (i in sct_lines) {
    if (!any(reset_lines > i)) {
      out[[length(out) + 1]] <- finding(
        "R-NORM-004", "L-NORM-002", "warn", i,
        'SCTransform() called with no later DefaultAssay(seu) <- "RNA" reset. DefaultAssay is now "SCT" - any later assay-implicit call silently operates on it.')
    }
  }

  if (!length(out)) return(no_findings())
  do.call(rbind, out)
}

format_findings <- function(f, source_name = "<input>") {
  if (!nrow(f)) return(sprintf("%s: clean\n", source_name))
  paste0(
    sprintf("%s: %d finding(s)\n", source_name, nrow(f)),
    paste(sprintf("  [%s] %s (%s) line %d\n    %s",
                  toupper(f$severity), f$rule, f$lint, f$line, f$message),
          collapse = "\n"), "\n")
}

if (!interactive() && identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) {
    worst <- 0L
    for (f in args) {
      res <- lint_norm_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
