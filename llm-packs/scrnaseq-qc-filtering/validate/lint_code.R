#!/usr/bin/env Rscript
#
# Static checks on R code produced by an LLM using the scrnaseq-qc-filtering
# pack. Deliberately dependency-free (no Seurat needed) so it can run in CI on a
# plain R image, and so an agent can self-check without a full environment.
#
# These are heuristics over source text, not semantic analysis. They are tuned
# for low false negatives on the failure modes in rules/qc-filtering.rules.yml;
# a few false positives are an acceptable price.
#
# Usage:
#   Rscript lint_code.R file1.R file2.R
#   lint_qc_code(text)   # when sourced

# --- helpers -----------------------------------------------------------------

strip_comments <- function(lines) {
  # Crude: drops anything after an unquoted #. Good enough - we only ever
  # over-strip, which loses context rather than inventing it.
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

# --- the checks --------------------------------------------------------------

lint_qc_code <- function(text) {
  lines <- if (length(text) == 1L) strsplit(text, "\n", fixed = TRUE)[[1]] else text
  code  <- strip_comments(lines)
  blob  <- paste(code, collapse = "\n")
  out   <- list()

  pfs_lines <- grep("PercentageFeatureSet", code)

  # L-001 (R-001): human mito pattern with no species determination nearby.
  #
  # Scoped to a window around the call, not the whole file. Whole-file scanning
  # produced a false negative on scripts/run_analysis.R, where the check was
  # satisfied by celldex::HumanPrimaryCellAtlasData() ~180 lines downstream -
  # a mention of "Human" with no bearing on the mitochondrial pattern.
  SPECIES_WINDOW <- 15L
  mt_upper <- grep('pattern\\s*=\\s*["\']\\^MT-', code)
  near_species <- function(i) {
    lo <- max(1L, i - SPECIES_WINDOW); hi <- min(length(code), i + SPECIES_WINDOW)
    any(grepl("species|organism|\\bmouse\\b|\\bhuman\\b|musculus|sapiens|melanogaster|switch\\s*\\(",
              code[lo:hi], ignore.case = TRUE))
  }
  species_aware <- length(mt_upper) && all(vapply(mt_upper, near_species, logical(1)))
  if (length(mt_upper) && !species_aware) {
    out[[length(out) + 1]] <- finding(
      "R-001", "L-001", "block", mt_upper[1],
      'Hardcoded "^MT-" with no species check. Matches zero genes on mouse data; percent.mito becomes 0 for every cell and filtering silently retains everything.')
  }

  # L-002 (R-002): PercentageFeatureSet without a non-empty-match guard.
  if (length(pfs_lines)) {
    guarded <- grepl("stopifnot|if\\s*\\(\\s*length\\s*\\(|length\\s*\\(\\s*grep", blob) &&
               grepl("grep\\s*\\(|rownames\\s*\\(", blob)
    if (!guarded) {
      out[[length(out) + 1]] <- finding(
        "R-002", "L-002", "block", pfs_lines[1],
        "PercentageFeatureSet called without asserting the pattern matched >0 features. An empty match returns 0 for all cells and does not warn.")
    }
  }

  # L-003 (R-004): Seurat 5 uses layer=, not slot=.
  slot_lines <- grep("slot\\s*=", code)
  if (length(slot_lines)) {
    out[[length(out) + 1]] <- finding(
      "R-004", "L-003", "warn", slot_lines[1],
      'slot= is deprecated in Seurat 5; use layer=. On a multi-sample v5 object also call JoinLayers() first, or you silently read one sample.')
  }

  # L-004 (R-003): over-broad haemoglobin pattern.
  hb_broad <- grep('pattern\\s*=\\s*["\']\\^HB["\']', code)
  if (length(hb_broad)) {
    out[[length(out) + 1]] <- finding(
      "R-003", "L-004", "warn", hb_broad[1],
      'Pattern "^HB" also matches HBEGF, HBS1L and HBP1, which are not haemoglobin subunits. Use an explicit gene list.')
  }

  # L-005 (R-010): a threshold LITERAL applied with no distribution inspected.
  #
  # Only numeric literals count. `percent.mito < max_mito` is a parameterised
  # function whose caller chooses the value and does the plotting - that is the
  # shape R-010 wants, not the shape it forbids. Matching it produced a false
  # positive on implementation/qc_filter.R:filter_cells().
  thresh_lines <- grep(
    "(percent\\.mito|percent\\.ribo|percent\\.globin|nFeature_RNA|nCount_RNA)\\s*[<>]=?\\s*[0-9]",
    code)
  plotted <- grepl("VlnPlot|FeatureScatter|ggplot|hist\\s*\\(|plot\\s*\\(", blob)
  if (length(thresh_lines) && !plotted) {
    out[[length(out) + 1]] <- finding(
      "R-010", "L-005", "warn", thresh_lines[1],
      "QC threshold applied without plotting the distribution first. Thresholds are dataset-specific decisions, not defaults.")
  }

  # L-006 (R-004): subset() after split() with no JoinLayers() in between.
  # Same failure family as L-003 (slot=/layer=), applied to subset() rather
  # than a GetAssayData()-style accessor - see reference/pitfalls.md#p-012.
  subset_lines <- grep("\\bsubset\\s*\\(", code)
  split_lines  <- grep("\\bsplit\\s*\\(", code)
  join_lines   <- grep("JoinLayers\\s*\\(", code)
  for (i in subset_lines) {
    if (any(split_lines < i) && !any(join_lines < i)) {
      out[[length(out) + 1]] <- finding(
        "R-004", "L-006", "warn", i,
        "subset() called after split() with no JoinLayers() in between. On an unjoined multi-sample v5 object this can silently filter within a single layer rather than consistently across every sample.")
    }
  }

  # L-007 (R-006): a named highly-expressed gene (e.g. MALAT1) removed with no
  # stated connection to the rest of the QC evidence nearby.
  MALAT_WINDOW <- 3L
  malat_lines <- grep("MALAT1", code)
  for (i in malat_lines) {
    looks_like_removal <- grepl(
      '!=\\s*["\']MALAT1["\']|setdiff\\s*\\(|%in%|exclude|drop', code[i], ignore.case = TRUE)
    if (looks_like_removal) {
      lo <- max(1L, i - MALAT_WINDOW); hi <- min(length(code), i + MALAT_WINDOW)
      # Justification is often written as a comment, so check the raw lines,
      # not the comment-stripped `code` used for pattern matching.
      justified <- grepl("stress|dying|QC|percent\\.mito|justif|because|evidence",
                        paste(lines[lo:hi], collapse = " "), ignore.case = TRUE)
      if (!justified) {
        out[[length(out) + 1]] <- finding(
          "R-006", "L-007", "warn", i,
          "MALAT1 appears to be removed with no stated connection to the rest of the QC evidence. State why (or leave it in, as the course chapter does when no other QC metric suggests stressed/dying cells).")
      }
    }
  }

  # L-008 (R-007): reading from the raw (not cellranger-filtered) matrix.
  raw_matrix_lines <- grep("raw_feature_bc_matrix", code)
  for (i in raw_matrix_lines) {
    out[[length(out) + 1]] <- finding(
      "R-007", "L-008", "warn", i,
      "Reading from raw_feature_bc_matrix/, not filtered_feature_bc_matrix/. This has not had cellranger's own empty-droplet cell-calling (the UMI-rank knee point) applied - it needs its own cell-calling step (e.g. emptyDrops()) before this pack's thresholds apply.")
  }

  # L-009 (R-009): PercentageFeatureSet() after NormalizeData()/SCTransform().
  norm_lines <- grep("NormalizeData\\s*\\(|SCTransform\\s*\\(", code)
  late_pfs <- Filter(function(i) any(norm_lines < i), pfs_lines)
  if (length(late_pfs)) {
    out[[length(out) + 1]] <- finding(
      "R-009", "L-009", "block", late_pfs[1],
      "PercentageFeatureSet() called after NormalizeData()/SCTransform(). QC metrics computed on normalised data are not a percentage of anything meaningful - compute them on raw counts, before normalisation.")
  }

  # L-010 (R-014): percent.ribo used as a filtering threshold, not just reported.
  ribo_filter <- grep("percent\\.ribo\\s*[<>]=?\\s*[0-9]", code)
  for (i in ribo_filter) {
    out[[length(out) + 1]] <- finding(
      "R-014", "L-010", "warn", i,
      "percent.ribo used as a filtering threshold. Unlike percent.mito, there is no widely-accepted mechanism linking ribosomal percentage to cell quality - report it, but don't filter on it without a stated reason for this dataset.")
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

# --- CLI ---------------------------------------------------------------------

if (!interactive() && identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) {
    worst <- 0L
    for (f in args) {
      res <- lint_qc_code(readLines(f, warn = FALSE))
      cat(format_findings(res, f))
      if (any(res$severity == "block")) worst <- 2L
      else if (nrow(res) && worst < 1L) worst <- 1L
    }
    quit(status = if (worst == 2L) 1L else 0L)
  }
}
