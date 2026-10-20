#!/usr/bin/env Rscript
#
# Static checks for scrnaseq-dimensionality-reduction.

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

lint_dimred_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  out   <- list()

  # L-DR-001 (R-DR-002): RunUMAP()/RunTSNE() with no evidence a PCA (or
  # PCA-derived reduction) exists first. Three kinds of evidence are
  # accepted, any one is enough:
  #   (a) an actual RunPCA() call earlier in the file
  #   (b) a runtime existence check against seu@reductions earlier (the
  #       reference implementation's own approach: check at call time rather
  #       than require a redundant RunPCA() in the same script)
  #   (c) a reduction= argument naming something other than a literal "pca"
  #       (e.g. an integrated reduction built on a PCA elsewhere)
  dimred_lines <- grep("RunUMAP\\s*\\(|RunTSNE\\s*\\(", code)
  pca_call_lines <- grep("\\bRunPCA\\s*\\(", code)
  reductions_mentioned <- grep("reductions", code)
  stop_lines <- grep("\\bstop\\s*\\(", code)
  for (i in dimred_lines) {
    names_other_reduction <- grepl('reduction\\s*=\\s*["\'](?!pca["\'])', code[i], perl = TRUE)
    has_prior_pca_call <- any(pca_call_lines < i)
    has_prior_runtime_check <- any(reductions_mentioned < i) && any(stop_lines < i)
    if (!has_prior_pca_call && !has_prior_runtime_check && !names_other_reduction) {
      out[[length(out) + 1]] <- finding(
        "R-DR-002", "L-DR-001", "block", i,
        "RunUMAP()/RunTSNE() called with no RunPCA() call, no runtime check that a PCA reduction exists, and no alternative reduction named. t-SNE/UMAP must not run before a PCA (or PCA-derived reduction) exists.")
    }
  }

  # L-DR-002 (R-DR-005): a literal dims = a:b range with no ElbowPlot (or
  # other PC-selection diagnostic) call anywhere in the file.
  dims_literal <- grep("dims\\s*=\\s*[0-9]+\\s*:\\s*[0-9]+", code)
  has_diagnostic <- any(grepl("ElbowPlot\\s*\\(|JackStraw", code))
  for (i in dims_literal) {
    if (!has_diagnostic) {
      out[[length(out) + 1]] <- finding(
        "R-DR-005", "L-DR-002", "warn", i,
        "A literal dims = a:b range with no ElbowPlot() (or other PC-selection diagnostic) anywhere in the file. State how this range was chosen for this object.")
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
      res <- lint_dimred_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
