#!/usr/bin/env Rscript
#
# Reference implementation: PCA-gated UMAP/t-SNE.
# Derived from day2/day2-2_dimensionality_reduction.qmd.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

#' Run UMAP (optionally t-SNE), refusing to proceed without a prior PCA
#'
#' @param seu       Seurat object with variable features set.
#' @param dims      PCs to use. No default - rule R-DR-005: this must follow
#'                  from an elbow plot (or equivalent) for THIS object.
#' @param reduction Which existing reduction to build on - "pca" by default,
#'                  or an integrated reduction name.
#' @param run_tsne  Also run t-SNE alongside UMAP.
run_dimred <- function(seu, dims, reduction = "pca", run_tsne = FALSE) {
  if (missing(dims)) {
    stop("`dims` is required (rule R-DR-005). Generate ElbowPlot(seu) for ",
         "THIS object and choose dims from where it flattens - a value ",
         "copied from another dataset is not a decision made for this one.")
  }
  # R-DR-002: PCA (or an equivalent prior reduction) must already exist.
  if (!reduction %in% names(seu@reductions)) {
    stop(sprintf(paste0(
      "Reduction '%s' not found on the object. Rule R-DR-002: t-SNE/UMAP ",
      "must not run on raw/scaled expression directly - run RunPCA() (or ",
      "the relevant integration step) first."), reduction))
  }
  available_pcs <- ncol(SeuratObject::Embeddings(seu, reduction))
  if (max(dims) > available_pcs) {
    stop(sprintf(paste0(
      "dims goes up to %d, but '%s' only has %d components. Increase npcs ",
      "at the RunPCA()/integration step if more precision is needed - ",
      "dims cannot exceed what was already computed."),
      max(dims), reduction, available_pcs))
  }

  seu <- Seurat::RunUMAP(seu, dims = dims, reduction = reduction, verbose = FALSE)
  if (run_tsne) {
    seu <- Seurat::RunTSNE(seu, dims = dims, reduction = reduction)
  }

  attr(seu, "dimred_dims")      <- dims
  attr(seu, "dimred_reduction") <- reduction
  attr(seu, "dimred_tsne")      <- run_tsne
  stamp_provenance(seu)
}

# --- Attribution -------------------------------------------------------------

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-dimensionality-reduction",
  version    = "0.1.0",
  title      = "Course material: Single cell transcriptomics with R",
  authors    = "Wyss T, Marcone R, van Geest G, Cascione L, Palagi P",
  publisher  = "SIB Swiss Institute of Bioinformatics",
  doi        = "10.5281/zenodo.5703106",
  url        = "https://github.com/sib-swiss/single-cell-r-training",
  license    = "CC-BY-4.0"
)

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Explicit arguments let this be called directly after an inline
#' RunPCA()/RunUMAP() sequence, using the same dims/reduction just passed to
#' them - not only after run_dimred().
stamp_provenance <- function(seu, dims = NULL, reduction = NULL, tsne = NULL) {
  dims_use <- dims %||% attr(seu, "dimred_dims")
  prov <- list(
    pack      = PACK_CITATION$pack,
    version   = PACK_CITATION$version,
    doi       = PACK_CITATION$doi,
    dims      = paste(range(dims_use), collapse = "-"),
    reduction = reduction %||% attr(seu, "dimred_reduction"),
    tsne      = if (!is.null(tsne)) tsne else attr(seu, "dimred_tsne"),
    run_at    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R         = as.character(getRversion()),
    Seurat    = as.character(utils::packageVersion("Seurat"))
  )
  existing <- attr(seu, "provenance")
  attr(seu, "provenance") <- c(existing, list(prov))
  seu
}

dimred_methods_sentence <- function(seu) {
  p <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  sprintf(paste0(
    "Dimensionality reduction was performed in Seurat v%s: UMAP%s was run on ",
    "dimensions %s of the '%s' reduction, following the %s workflow (%s; ",
    "%s, %s, DOI: %s)."),
    p$Seurat, if (isTRUE(p$tsne)) " and t-SNE" else "", p$dims, p$reduction,
    PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
    PACK_CITATION$publisher, p$doi)
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
