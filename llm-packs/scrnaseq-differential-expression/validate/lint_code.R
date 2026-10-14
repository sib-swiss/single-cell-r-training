#!/usr/bin/env Rscript
#
# Static checks on R code produced by an LLM using the
# scrnaseq-differential-expression pack. Dependency-free (no Seurat needed),
# heuristic over source text - same philosophy as the QC pack's linter.
#
# Usage:
#   Rscript lint_code.R file1.R file2.R
#   lint_de_code(text)   # when sourced

CONDITION_KEYWORDS <- c("condition", "type", "group", "treatment", "status",
                        "disease", "genotype", "stim", "stimulated")
REPLICATE_KEYWORDS <- c("sample", "orig.ident", "patient", "donor", "subject",
                        "mouse", "animal", "individual")

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

keyword_pattern <- function(words) paste0("(", paste(words, collapse = "|"), ")")

lint_de_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  blob  <- paste(code, collapse = "\n")
  out   <- list()

  has_pseudobulk <- grepl("AggregateExpression\\s*\\(|pseudobulk_de\\s*\\(", blob, ignore.case = TRUE)

  # L-DE-001 (R-DE-001): a cross-condition FindMarkers/FindMarkers-style call
  # (group.by/ident args reference a condition-like keyword) with no
  # pseudobulk aggregation anywhere in the file.
  de_lines <- grep("FindMarkers\\s*\\(|FindAllMarkers\\s*\\(", code)
  for (i in de_lines) {
    window <- paste(code[max(1, i - 3):min(length(code), i + 3)], collapse = " ")
    looks_cross_condition <- grepl(
      paste0("group\\.by\\s*=\\s*[^,\\)]*", keyword_pattern(CONDITION_KEYWORDS)),
      window, ignore.case = TRUE)
    if (looks_cross_condition && !has_pseudobulk) {
      out[[length(out) + 1]] <- finding(
        "R-DE-001", "L-DE-001", "block", i,
        "FindMarkers/FindAllMarkers appears to compare a condition that spans multiple samples, with no pseudobulk aggregation (AggregateExpression/pseudobulk_de) anywhere in the file. Per-cell tests across patients/samples suffer pseudoreplication - route cross-condition comparisons through pseudobulk aggregation first.")
    }
  }

  # L-DE-002 (R-DE-003): AggregateExpression() grouped only by a
  # condition-like variable, with no replicate/sample-like keyword in the
  # group.by value - pools multiple replicates of the same condition into one
  # pseudobulk profile, defeating the purpose of pseudobulk testing.
  agg_lines <- grep("AggregateExpression\\s*\\(", code)
  for (i in agg_lines) {
    window <- paste(code[i:min(length(code), i + 4)], collapse = " ")
    gb <- regmatches(window, regexec('group\\.by\\s*=\\s*["\']([^"\']+)["\']', window))[[1]]
    if (length(gb) == 2) {
      val <- gb[2]
      is_condition_only <- grepl(keyword_pattern(CONDITION_KEYWORDS), val, ignore.case = TRUE) &&
                           !grepl(keyword_pattern(REPLICATE_KEYWORDS), val, ignore.case = TRUE)
      if (is_condition_only) {
        out[[length(out) + 1]] <- finding(
          "R-DE-003", "L-DE-002", "warn", i,
          sprintf("AggregateExpression(group.by = \"%s\") groups only by condition, not by sample/replicate. This pools every biological replicate of a condition into one pseudobulk profile, which removes the ability to estimate between-replicate variance - group by the sample/patient identifier instead.", val))
      }
    }
  }

  # L-DE-003 (R-DE-005): a DE result is filtered/reported as significant with
  # no diagnostic plot (p-value histogram or MA-style plot) anywhere in the
  # file.
  sig_filter <- grep("p_val_adj\\s*<|adj\\.P\\.Val\\s*<", code)
  has_diagnostic <- grepl("geom_histogram|hist\\s*\\(|de_diagnostic_plots\\s*\\(", blob)
  if (length(sig_filter) && !has_diagnostic) {
    out[[length(out) + 1]] <- finding(
      "R-DE-005", "L-DE-003", "warn", sig_filter[1],
      "Genes filtered to 'significant' with no diagnostic plot (p-value histogram or MA plot / de_diagnostic_plots()) anywhere in the file. Inspect these before reporting the gene list as trustworthy.")
  }

  # L-DE-005 (R-DE-004): 2+ separate pairwise FindMarkers calls (implying
  # >2 groups compared pairwise) with no omnibus/ANOVA alternative in sight.
  if (length(de_lines) >= 2) {
    has_omnibus <- grepl("anova_de\\s*\\(|\\bF-test\\b|topTable\\s*\\([^)]*\\)\\s*(?!.*coef)",
                        blob, ignore.case = TRUE, perl = TRUE)
    if (!has_omnibus) {
      out[[length(out) + 1]] <- finding(
        "R-DE-004", "L-DE-005", "warn", de_lines[1],
        "Multiple separate FindMarkers calls found (pairwise comparisons across what may be >2 groups), with no omnibus/ANOVA-style test (e.g. anova_de(), a limma F-test) in the file. Each pairwise call's p_val_adj is only adjusted within that comparison - this under-corrects for testing across all pairs. Consider an omnibus test, or explicitly re-adjust the combined p-values.")
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
      res <- lint_de_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
