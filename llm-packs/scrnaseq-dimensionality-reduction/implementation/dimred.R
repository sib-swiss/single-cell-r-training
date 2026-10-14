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
