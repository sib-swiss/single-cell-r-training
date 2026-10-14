#!/usr/bin/env Rscript
#
# Static checks on R code produced by an LLM using the scrnaseq-integration
# pack. Dependency-free (no Seurat needed), heuristic over source text rather
# than semantic analysis - same philosophy as the QC pack's linter.
#
# Usage:
#   Rscript lint_code.R file1.R file2.R
#   lint_integration_code(text)   # when sourced

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

lint_integration_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  out   <- list()

  # L-INT-001 (R-INT-001): DE/marker calls on the legacy "integrated" assay.
  #
  # Track the default assay as we scan sequentially - the same tracked-state
  # approach is needed here because the dangerous pattern spans two separate
  # statements (DefaultAssay<- then FindMarkers), unlike the QC pack's
  # single-line findings.
  current_assay <- NA_character_
  de_calls <- grep("FindMarkers|FindAllMarkers", code)
  for (i in seq_along(code)) {
    m <- regmatches(code[i], regexec('DefaultAssay\\s*\\([^)]*\\)\\s*<-\\s*["\']([^"\']+)["\']', code[i]))[[1]]
    if (length(m) == 2) current_assay <- m[2]
    if (i %in% de_calls) {
      explicit_assay <- regmatches(code[i], regexec('assay\\s*=\\s*["\']([^"\']+)["\']', code[i]))[[1]]
      effective_assay <- if (length(explicit_assay) == 2) explicit_assay[2] else current_assay
      if (!is.na(effective_assay) && tolower(effective_assay) == "integrated") {
        out[[length(out) + 1]] <- finding(
          "R-INT-001", "L-INT-001", "block", i,
          'FindMarkers/FindAllMarkers called with the default assay set to "integrated". IntegrateLayers/IntegrateData output is corrected for visualisation and clustering, not for DE - it gives the batch more statistical weight than the biology. Switch back to the original (joined, non-integrated) assay before testing.')
      }
    }
  }

  # Legacy anchor-based workflow (IntegrateData) is the one that actually
  # produces a DE-unsafe "integrated" assay by default - flag its use outright
  # so the model considers the v5 layer-based workflow instead.
  legacy <- grep("\\bIntegrateData\\s*\\(", code)
  if (length(legacy)) {
    out[[length(out) + 1]] <- finding(
      "R-INT-001", "L-INT-001", "warn", legacy[1],
      "IntegrateData() is the pre-v5 anchor-based workflow and produces a corrected 'integrated' assay intended for clustering/visualisation only. This course teaches the v5 layer-based IntegrateLayers() workflow instead; if IntegrateData() is intentional, the resulting assay must still not be used for DE (R-INT-001).")
  }

  # L-INT-002 (R-INT-004): IntegrateLayers() with no JoinLayers() anywhere
  # after it in the file.
  int_lines <- grep("IntegrateLayers\\s*\\(", code)
  join_lines <- grep("JoinLayers\\s*\\(", code)
  for (i in int_lines) {
    if (!any(join_lines > i)) {
      out[[length(out) + 1]] <- finding(
        "R-INT-004", "L-INT-002", "warn", i,
        "IntegrateLayers() called with no later JoinLayers() call in this file. A multi-sample v5 object keeps per-sample layers until joined; any single-matrix read before that silently covers one sample.")
    }
  }

  # L-INT-003 (R-INT-005): SCT-assay integration reusing the RNA assay's
  # plain "pca" reduction. Scoped to the IntegrateLayers(...) call text
  # (this call, plus the following few lines, to tolerate multi-line calls)
  # rather than the whole file, for the same reason the QC pack's L-001
  # scopes to a window: a distant, unrelated mention must not suppress or
  # trigger a finding.
  CALL_WINDOW <- 8L
  for (i in int_lines) {
    hi <- min(length(code), i + CALL_WINDOW)
    call_text <- paste(code[i:hi], collapse = " ")
    if (grepl('assay\\s*=\\s*["\']SCT', call_text) &&
        grepl('orig\\.reduction\\s*=\\s*["\']pca["\']', call_text)) {
      out[[length(out) + 1]] <- finding(
        "R-INT-005", "L-INT-003", "block", i,
        'IntegrateLayers() integrating the SCT assay with orig.reduction = "pca" (the RNA assay\'s PCA). Compute a PCA on the SCT assay itself (e.g. reduction.name = "pcaSCT") and reference that instead - reusing the RNA PCA integrates on the wrong basis without erroring.')
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
      res <- lint_integration_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
