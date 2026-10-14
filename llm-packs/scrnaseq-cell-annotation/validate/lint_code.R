#!/usr/bin/env Rscript
#
# Static checks for scrnaseq-cell-annotation.

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

lint_annotation_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  blob  <- paste(code, collapse = "\n")
  out   <- list()

  # L-ANN-001 (R-ANN-001): slot= is deprecated in Seurat 5.
  slot_lines <- grep("slot\\s*=", code)
  for (i in slot_lines) {
    out[[length(out) + 1]] <- finding(
      "R-ANN-001", "L-ANN-001", "warn", i,
      'slot= is deprecated in Seurat 5; use layer= (or omit it - GetAssayData() defaults to layer = "data"). Known live instance: scripts/run_analysis.R uses slot= for the SingleR test matrix.')
  }

  # L-ANN-002 (R-ANN-002): SingleR() called with no diagnostic plot anywhere.
  singler_lines <- grep("SingleR::SingleR\\s*\\(|\\bSingleR\\s*\\(", code)
  has_diagnostic <- grepl("plotScoreHeatmap\\s*\\(|plotDeltaDistribution\\s*\\(|annotation_diagnostic_plots\\s*\\(", blob)
  if (length(singler_lines) && !has_diagnostic) {
    out[[length(out) + 1]] <- finding(
      "R-ANN-002", "L-ANN-002", "warn", singler_lines[1],
      "SingleR() called with no plotScoreHeatmap()/plotDeltaDistribution() diagnostic anywhere in the file. SingleR has no rejection option - every cell gets a label - inspect quality before trusting the annotation.")
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
      res <- lint_annotation_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
