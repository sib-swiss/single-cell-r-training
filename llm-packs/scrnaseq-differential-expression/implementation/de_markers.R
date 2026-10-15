#!/usr/bin/env Rscript
#
# Reference implementation: differential gene expression, Seurat 5 + edgeR/limma.
#
# This is the artifact an LLM should ADAPT rather than reimplement. It is
# exercised by validate/run_all.R in CI against synthetic fixtures, so its
# behaviour on the failure modes in rules/differential-expression.rules.yml is
# a tested property, not a claim.
#
# Derived from day3/day3-2_differential_gene_expression.qmd, with the guards
# the teaching version omits for readability.

suppressPackageStartupMessages({
  library(Seurat)
})

# --- Marker discovery (per-cluster, exploratory) -----------------------------
#
# Per-cell testing is appropriate here (R-DE-001 only concerns CROSS-CONDITION
# comparisons spanning multiple samples): cluster membership is a per-cell
# property, not a per-patient one, so there is no pseudoreplication problem in
# the same sense.

#' Find marker genes for every cluster (one cluster vs. the rest)
cluster_markers <- function(seu, min_pct = 0.25, only_pos = TRUE,
                            test.use = "wilcox") {
  markers <- Seurat::FindAllMarkers(seu, min.pct = min_pct, only.pos = only_pos,
                                    test.use = test.use, verbose = FALSE)
  markers <- subset(markers, markers$p_val_adj < 0.05)
  stamp_provenance(markers, step = "cluster_markers",
                   n_clusters = length(unique(markers$cluster)))
}

# --- Cross-condition comparison (pseudobulk) ---------------------------------

#' Aggregate single cells to one pseudobulk profile per (replicate x condition)
#' and test for differential expression with edgeR + limma-voom.
#'
#' @param seu            Seurat object.
#' @param sample_col     metadata column identifying the sample (the unit that
#'                       gets aggregated into one pseudobulk column).
#' @param condition_col  metadata column giving the condition to compare.
#' @param replicate_col  metadata column identifying the TRUE independent
#'                       biological replicate (e.g. patient/animal ID), if a
#'                       paired/blocked design is intended. NULL for an
#'                       unpaired design. No default - forcing the caller to
#'                       decide is rule R-DE-002: there is no way to verify
#'                       from the data alone whether a replicate_col reflects
#'                       real pairing, so the caller must have confirmed it.
#' @param contrast       length-2 character vector: c(level_of_interest, reference_level).
#' @param assay          Assay to aggregate from.
#' @param min_replicates Minimum independent replicates required per
#'                       condition (rule R-DE-001 / pack.yml:not_applicable_to).
pseudobulk_de <- function(seu, sample_col, condition_col, contrast,
                          replicate_col = NULL, assay = "RNA",
                          min_replicates = 2) {
  if (missing(sample_col) || missing(condition_col) || missing(contrast)) {
    stop("`sample_col`, `condition_col` and `contrast` are all required.")
  }
  md <- seu@meta.data
  n_per_condition <- table(unique(md[, c(sample_col, condition_col)])[[condition_col]])
  if (any(n_per_condition[contrast] < min_replicates)) {
    stop(sprintf(paste0(
      "Fewer than %d independent replicates (distinct '%s' values) in at ",
      "least one of %s. Rule R-DE-001: a cross-condition comparison with ",
      "this few replicates cannot distinguish biological variability from ",
      "noise, and must not be run as a per-cell test to compensate - there ",
      "is no substitute for more replicates."),
      min_replicates, sample_col, paste(contrast, collapse = " vs ")))
  }

  if (!is.null(replicate_col)) {
    # Structural sanity check only - this CANNOT verify real-world pairing
    # (rule R-DE-002 is a manual check), it can only catch an inconsistent
    # mapping.
    map <- unique(md[, c(sample_col, replicate_col, condition_col)])
    dup <- duplicated(map[[sample_col]])
    if (any(dup)) {
      stop("Each value of `sample_col` must map to exactly one `replicate_col` ",
           "and one `condition_col` - found an inconsistent mapping.")
    }
    message("`replicate_col` = '", replicate_col, "' implies a paired design. ",
            "This has NOT been verified against real experimental metadata - ",
            "confirm with the user that samples are genuinely matched before ",
            "reporting a paired result (rule R-DE-002).")
  }

  bulk <- Seurat::AggregateExpression(seu, group.by = sample_col,
                                      return.seurat = TRUE, assay = assay,
                                      verbose = FALSE)
  meta <- unique(md[, c(sample_col, condition_col, replicate_col)])
  rownames(meta) <- meta[[sample_col]]
  bulk@meta.data <- meta[colnames(bulk), , drop = FALSE]

  counts <- as.matrix(Seurat::GetAssayData(bulk, layer = "counts"))
  y <- edgeR::DGEList(counts, samples = bulk@meta.data)
  keep <- edgeR::filterByExpr(y, group = bulk@meta.data[[condition_col]])
  y <- y[keep, ]

  dat <- bulk@meta.data
  dat$.condition <- factor(dat[[condition_col]])
  if (!is.null(replicate_col)) dat$.replicate <- factor(dat[[replicate_col]])

  form <- if (is.null(replicate_col)) ~0 + .condition else ~0 + .condition + .replicate
  design <- stats::model.matrix(form, data = dat)
  colnames(design) <- sub("^\\.condition", "", colnames(design))
  colnames(design) <- sub("^\\.replicate", "rep", colnames(design))
  colnames(design) <- make.names(colnames(design))

  lvl_names <- make.names(contrast)
  contrast_str <- sprintf("%s - %s", lvl_names[1], lvl_names[2])
  cm <- limma::makeContrasts(contrasts = contrast_str, levels = design)

  v <- limma::voom(edgeR::calcNormFactors(y), design)
  fit <- limma::contrasts.fit(limma::lmFit(v, design), cm)
  fit <- limma::eBayes(fit)

  res <- limma::topTable(fit, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res <- res[, c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")]

  stamp_provenance(res, step = "pseudobulk_de",
                   contrast = paste(contrast, collapse = " vs "),
                   n_replicates = as.list(n_per_condition[contrast]),
                   paired = !is.null(replicate_col))
}

#' Omnibus (ANOVA-style) test across 3+ groups, correctly controlling the
#' family-wise error rate - the R-DE-004 alternative to repeated pairwise
#' FindMarkers calls.
anova_de <- function(seu, sample_col, group_col, assay = "RNA") {
  md <- seu@meta.data
  groups <- unique(md[[group_col]])
  if (length(groups) < 3) {
    stop("anova_de() is for 3+ groups. For 2 groups use pseudobulk_de().")
  }
  bulk <- Seurat::AggregateExpression(seu, group.by = sample_col,
                                      return.seurat = TRUE, assay = assay,
                                      verbose = FALSE)
  meta <- unique(md[, c(sample_col, group_col)])
  rownames(meta) <- meta[[sample_col]]
  bulk@meta.data <- meta[colnames(bulk), , drop = FALSE]

  counts <- as.matrix(Seurat::GetAssayData(bulk, layer = "counts"))
  y <- edgeR::DGEList(counts, samples = bulk@meta.data)
  y <- y[edgeR::filterByExpr(y, group = bulk@meta.data[[group_col]]), ]

  design <- stats::model.matrix(stats::as.formula(sprintf("~0 + %s", group_col)),
                                data = bulk@meta.data)
  colnames(design) <- make.names(sub(paste0("^", group_col), "", colnames(design)))

  v <- limma::voom(edgeR::calcNormFactors(y), design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  # F-test across all coefficients at once - this IS the correctly-adjusted
  # omnibus test rule R-DE-004 asks for, in place of repeated pairwise calls.
  res <- limma::topTable(fit, number = Inf, sort.by = "F")
  res$gene <- rownames(res)
  stamp_provenance(res, step = "anova_de", groups = paste(groups, collapse = ", "))
}

# --- Diagnostics (R-DE-005) ---------------------------------------------------

#' p-value histogram + MA plot. A DE result should not be reported as
#' trustworthy without having looked at both.
de_diagnostic_plots <- function(de_result) {
  pcol <- if ("adj.P.Val" %in% names(de_result)) "adj.P.Val" else "p_val_adj"
  fccol <- if ("logFC" %in% names(de_result)) "logFC" else "avg_log2FC"
  avcol <- if ("AveExpr" %in% names(de_result)) "AveExpr" else NA

  p_hist <- ggplot2::ggplot(de_result, ggplot2::aes(.data[[pcol]])) +
    ggplot2::geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "p-value histogram (expect roughly L-shaped)",
                 x = pcol, y = "count") +
    ggplot2::theme_minimal()

  ma <- if (!is.na(avcol)) {
    ggplot2::ggplot(de_result, ggplot2::aes(.data[[avcol]], .data[[fccol]])) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::geom_hline(yintercept = 0, color = "black") +
      ggplot2::labs(title = "MA plot (expect symmetry around logFC = 0)",
                   x = avcol, y = fccol) +
      ggplot2::theme_minimal()
  } else NULL

  list(p_value_histogram = p_hist, ma_plot = ma)
}

# --- Attribution -------------------------------------------------------------

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-differential-expression",
  version    = "0.1.0",
  title      = "Course material: Single cell transcriptomics with R",
  authors    = "Wyss T, Marcone R, van Geest G, Cascione L, Palagi P",
  publisher  = "SIB Swiss Institute of Bioinformatics",
  doi        = "10.5281/zenodo.5703106",
  url        = "https://github.com/sib-swiss/single-cell-r-training",
  license    = "CC-BY-4.0"
)

#' Attach a provenance record as an attribute on a data.frame result
stamp_provenance <- function(df, step, ...) {
  prov <- c(list(
    pack    = PACK_CITATION$pack,
    version = PACK_CITATION$version,
    doi     = PACK_CITATION$doi,
    step    = step,
    run_at  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R       = as.character(getRversion())
  ), list(...))
  existing <- attr(df, "provenance")
  attr(df, "provenance") <- c(existing, list(prov))
  df
}

#' A methods sentence the researcher can paste into a manuscript
de_methods_sentence <- function(de_result) {
  p <- utils::tail(attr(de_result, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  if (identical(p$step, "pseudobulk_de")) {
    sprintf(paste0(
      "Differential expression between %s was assessed on pseudobulk ",
      "profiles (one per independent biological replicate) using edgeR/",
      "limma-voom, following the %s workflow (%s; %s, %s, DOI: %s). %s"),
      p$contrast, PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
      PACK_CITATION$publisher, p$doi,
      if (isTRUE(p$paired)) "A paired design was used." else "An unpaired design was used.")
  } else {
    sprintf(paste0(
      "Marker genes were identified using Seurat's FindAllMarkers (Wilcoxon ",
      "test) following the %s workflow (%s; %s, %s, DOI: %s)."),
      PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
      PACK_CITATION$publisher, p$doi)
  }
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

write_provenance <- function(de_result, path = "provenance.json") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not available; provenance not written"); return(invisible(NULL))
  }
  jsonlite::write_json(list(provenance = attr(de_result, "provenance"),
                            citation = PACK_CITATION),
                       path, auto_unbox = TRUE, pretty = TRUE)
  message("Provenance written to ", path)
  invisible(path)
}
