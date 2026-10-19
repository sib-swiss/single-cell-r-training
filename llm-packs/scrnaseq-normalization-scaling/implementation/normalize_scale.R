#!/usr/bin/env Rscript
#
# Reference implementation: normalization, variable features and scaling.
# Derived from day2/day2-1_normalization_scaling.qmd.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

#' Normalize, select variable features, and scale
#'
#' @param seu              QC-filtered Seurat object.
#' @param nfeatures        Variable feature count. Has a default (2000, the
#'                         course's own convention) unlike the QC pack's
#'                         thresholds, because this is a computational
#'                         trade-off knob, not a biological decision
#'                         pretending to be universal - but R-NORM-002 still
#'                         requires stating it as a convention, not fact.
#' @param scale_factor     NormalizeData's scale.factor. Convention: must be
#'                         large relative to typical per-cell library size.
#' @param extra_scale_genes Genes to scale in ADDITION to the variable
#'                         features (R-NORM-001) - e.g. a marker panel a
#'                         downstream heatmap will need.
scale_variable_features <- function(seu, nfeatures = 2000, scale_factor = 10000,
                                    extra_scale_genes = character(0)) {
  seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                               scale.factor = scale_factor, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                      nfeatures = nfeatures, verbose = FALSE)

  scale_genes <- union(Seurat::VariableFeatures(seu),
                       intersect(extra_scale_genes, rownames(seu)))
  missing_extra <- setdiff(extra_scale_genes, rownames(seu))
  if (length(missing_extra)) {
    warning("extra_scale_genes not found in this object and skipped: ",
            paste(missing_extra, collapse = ", "))
  }
  seu <- Seurat::ScaleData(seu, features = scale_genes, verbose = FALSE)

  attr(seu, "norm_nfeatures")    <- nfeatures
  attr(seu, "norm_scale_factor") <- scale_factor
  attr(seu, "norm_scaled_genes") <- length(scale_genes)
  stamp_provenance(seu)
}

#' Switch to SCTransform, warning about the DefaultAssay change (R-NORM-004)
run_sctransform <- function(seu, ...) {
  seu <- Seurat::SCTransform(seu, ...)
  message("DefaultAssay(seu) is now 'SCT'. Reset to 'RNA' explicitly ",
          "(DefaultAssay(seu) <- 'RNA') before any step that expects the ",
          "RNA assay - rule R-NORM-004.")
  attr(seu, "norm_method") <- "SCTransform"
  stamp_provenance(seu)
}

# --- Attribution -------------------------------------------------------------

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-normalization-scaling",
  version    = "0.1.0",
  title      = "Course material: Single cell transcriptomics with R",
  authors    = "Wyss T, Marcone R, van Geest G, Cascione L, Palagi P",
  publisher  = "SIB Swiss Institute of Bioinformatics",
  doi        = "10.5281/zenodo.5703106",
  url        = "https://github.com/sib-swiss/single-cell-r-training",
  license    = "CC-BY-4.0"
)

#' Explicit arguments let this be called directly after inline
#' NormalizeData()/FindVariableFeatures()/ScaleData() calls, using the same
#' values just passed to them - not only after scale_variable_features().
stamp_provenance <- function(seu, method = NULL, nfeatures = NULL,
                             scale_factor = NULL, n_scaled_genes = NULL) {
  prov <- list(
    pack        = PACK_CITATION$pack,
    version     = PACK_CITATION$version,
    doi         = PACK_CITATION$doi,
    method      = method %||% attr(seu, "norm_method") %||% "LogNormalize",
    nfeatures   = nfeatures %||% attr(seu, "norm_nfeatures"),
    scale_factor = scale_factor %||% attr(seu, "norm_scale_factor"),
    n_scaled_genes = n_scaled_genes %||% attr(seu, "norm_scaled_genes"),
    run_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R           = as.character(getRversion()),
    Seurat      = as.character(utils::packageVersion("Seurat"))
  )
  existing <- attr(seu, "provenance")
  attr(seu, "provenance") <- c(existing, list(prov))
  seu
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

normalization_methods_sentence <- function(seu) {
  p <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  if (identical(p$method, "SCTransform")) {
    sprintf(paste0(
      "Data were normalised with SCTransform in Seurat v%s, following the %s ",
      "workflow (%s; %s, %s, DOI: %s)."),
      p$Seurat, PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
      PACK_CITATION$publisher, p$doi)
  } else {
    sprintf(paste0(
      "Data were log-normalised (scale factor %d) in Seurat v%s, with %d ",
      "variable features selected and %d genes scaled, following the %s ",
      "workflow (%s; %s, %s, DOI: %s)."),
      p$scale_factor, p$Seurat, p$nfeatures, p$n_scaled_genes,
      PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
      PACK_CITATION$publisher, p$doi)
  }
}

write_provenance <- function(seu, path = "provenance.json") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not available; provenance not written"); return(invisible(NULL))
  }
  jsonlite::write_json(list(provenance = attr(seu, "provenance"),
                            citation = PACK_CITATION),
                       path, auto_unbox = TRUE, pretty = TRUE)
  message("Provenance written to ", path)
  invisible(path)
}

#' Append (or update) this step's entry in the project's running methods
#' file - the cumulative, human-readable record across every pack used in
#' this project, maintained as the analysis moves forward (provenance.json
#' is the per-object machine record; this is the manuscript-facing one).
#' Keyed by `section`, so re-running the same step replaces its own entry
#' instead of accumulating stale duplicates.
#'
#' @param sentence A methods sentence, e.g. from *_methods_sentence().
#' @param path     Methods file path. Defaults to here::here("METHODS.md")
#'                 when the `here` package is available, else "METHODS.md"
#'                 in the working directory.
#' @param section  Key identifying this entry - a pack id, or
#'                 "<pack id>:<variant>" for a pack with more than one
#'                 workflow (e.g. differential expression's cluster_markers
#'                 vs pseudobulk_de).
append_methods <- function(sentence, path = NULL, section = PACK_CITATION$pack) {
  if (is.null(path)) {
    path <- if (requireNamespace("here", quietly = TRUE)) here::here("METHODS.md") else "METHODS.md"
  }
  begin_marker <- sprintf("<!-- methods:%s:begin -->", section)
  end_marker   <- sprintf("<!-- methods:%s:end -->", section)
  block <- c(begin_marker, "", sentence, "", end_marker)

  if (!file.exists(path)) {
    writeLines(c(
      "# Methods", "",
      "Built up automatically as each LLM pack step ran in this project.",
      "Each section below is replaced, not duplicated, the next time that",
      "step is re-run. Review before use - see each pack's",
      "reference/citations.md for what is independently verified versus",
      "sourced to course material.", ""), path)
  }

  lines <- readLines(path, warn = FALSE)
  start_idx <- which(lines == begin_marker)
  end_idx   <- which(lines == end_marker)
  if (length(start_idx) && length(end_idx)) {
    lines <- c(lines[seq_len(start_idx[1] - 1)], block,
              lines[seq(end_idx[1] + 1, length(lines))])
  } else {
    lines <- c(lines, "", block)
  }
  writeLines(lines, path)
  invisible(path)
}
