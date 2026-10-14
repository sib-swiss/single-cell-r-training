#!/usr/bin/env Rscript
#
# Static checks for scrnaseq-enrichment-analysis.

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

lint_enrichment_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  blob  <- paste(code, collapse = "\n")
  out   <- list()

  # L-ENR-001 (R-ENR-001): ORA call with no universe= argument. Scoped to a
  # small window to tolerate a multi-line call.
  ora_lines <- grep("enrichGO\\s*\\(|enricher\\s*\\(|enrichKEGG\\s*\\(", code)
  CALL_WINDOW <- 10L
  for (i in ora_lines) {
    hi <- min(length(code), i + CALL_WINDOW)
    call_text <- paste(code[i:hi], collapse = " ")
    if (!grepl("universe\\s*=", call_text)) {
      out[[length(out) + 1]] <- finding(
        "R-ENR-001", "L-ENR-001", "block", i,
        "ORA call with no universe= argument. Without it, clusterProfiler defaults to the whole annotated genome as background, not the genes actually tested in this experiment.")
    }
  }

  # L-ENR-002 (R-ENR-003): a GSEA-family call with no sort()/order() call
  # earlier in the file - a proxy for "was this list actually ranked first".
  gsea_lines <- grep("\\bGSEA\\s*\\(|gseGO\\s*\\(|gseKEGG\\s*\\(", code)
  sort_lines <- grep("\\bsort\\s*\\(|\\border\\s*\\(", code)
  for (i in gsea_lines) {
    if (!any(sort_lines < i)) {
      out[[length(out) + 1]] <- finding(
        "R-ENR-003", "L-ENR-002", "block", i,
        "GSEA-family call with no sort()/order() call earlier in the file. GSEA requires the input pre-ranked by a statistic - it does not rank the gene list for you.")
    }
  }

  # L-ENR-003 (R-ENR-005): an enrichGO() result's @result accessed/reported
  # with no simplify() call anywhere in the file.
  go_lines <- grep("enrichGO\\s*\\(", code)
  result_access <- grep("@result\\b", code)
  has_simplify <- grepl("\\bsimplify\\s*\\(", blob)
  if (length(go_lines) && length(result_access) && !has_simplify) {
    out[[length(out) + 1]] <- finding(
      "R-ENR-005", "L-ENR-003", "inform", result_access[1],
      "enrichGO() result accessed with no simplify() call anywhere in the file. GO terms are highly redundant - consider clusterProfiler::simplify() before reporting results as independent findings.")
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
      res <- lint_enrichment_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
