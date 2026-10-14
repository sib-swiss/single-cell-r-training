#!/usr/bin/env Rscript
#
# Static checks for scrnaseq-clustering.

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

lint_clustering_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  blob  <- paste(code, collapse = "\n")
  out   <- list()

  # L-CLU-001 (R-CLU-001): FindClusters() with a single literal resolution
  # (not a sweep/vector) and no clustree() stability check anywhere.
  single_res <- grep("FindClusters\\s*\\([^)]*resolution\\s*=\\s*[0-9.]+\\s*[,)]", code)
  has_clustree <- grepl("clustree\\s*\\(", blob)
  for (i in single_res) {
    if (!has_clustree) {
      out[[length(out) + 1]] <- finding(
        "R-CLU-001", "L-CLU-001", "warn", i,
        "FindClusters() called with a single literal resolution and no clustree() stability check anywhere in the file. Sweep a range of resolutions and check stability before presenting one as the result.")
    }
  }

  # L-CLU-002 (R-CLU-002): FindNeighbors() with no reduction= argument, in a
  # file that also shows evidence of integration.
  fn_lines <- grep("FindNeighbors\\s*\\(", code)
  looks_integrated <- grepl("IntegrateLayers\\s*\\(|integrated\\.cca|[Hh]armony|IntegrateData\\s*\\(",
                            blob)
  for (i in fn_lines) {
    has_reduction_arg <- grepl("reduction\\s*=", code[i])
    if (looks_integrated && !has_reduction_arg) {
      out[[length(out) + 1]] <- finding(
        "R-CLU-002", "L-CLU-002", "block", i,
        "FindNeighbors() called with no reduction= argument, in a file that shows evidence of integration. The default ('pca') would cluster on the unintegrated embedding - set reduction= to the integrated reduction's name explicitly.")
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
      res <- lint_clustering_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
