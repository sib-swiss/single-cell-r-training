#!/usr/bin/env Rscript
#
# Reference implementation: SNN graph + Louvain clustering.
# Derived from day2/day2-4_clustering.qmd.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

#' Build the neighbour graph and cluster across a resolution sweep
#'
#' @param seu         Seurat object with a PCA (and, if multi-sample, an
#'                    integrated) reduction.
#' @param dims        PCs to use.
#' @param reduction   Which reduction to build the graph on. No default when
#'                    an integrated/harmony/cca reduction is present on the
#'                    object - rule R-CLU-002: FindNeighbors()'s own default
#'                    ("pca") would silently ignore the batch correction.
#' @param resolutions Resolutions to sweep. Deliberately a sweep, not one
#'                    value - rule R-CLU-001: picking a single resolution
#'                    with no stability check is the thing this pack warns
#'                    against.
cluster_cells <- function(seu, dims, reduction = NULL,
                          resolutions = seq(0.1, 0.8, by = 0.1)) {
  if (missing(dims)) stop("`dims` is required.")

  integrated_reductions <- grep("integrated|harmony|cca|mnn", names(seu@reductions),
                                ignore.case = TRUE, value = TRUE)
  if (is.null(reduction)) {
    if (length(integrated_reductions)) {
      stop(sprintf(paste0(
        "This object carries an integrated reduction (%s), but `reduction` ",
        "was not specified. Rule R-CLU-002: the neighbour-graph step's own ",
        "default ('pca') would silently cluster on the UNINTEGRATED embedding. Pass ",
        "reduction = '%s' (or whichever integrated reduction is correct) ",
        "explicitly."), paste(integrated_reductions, collapse = ", "),
        integrated_reductions[1]))
    }
    reduction <- "pca"
  }
  if (!reduction %in% names(seu@reductions)) {
    stop(sprintf("Reduction '%s' not found on the object.", reduction))
  }

  seu <- Seurat::FindNeighbors(seu, dims = dims, reduction = reduction, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = resolutions, verbose = FALSE)

  attr(seu, "clustering_reduction")  <- reduction
  attr(seu, "clustering_resolutions") <- resolutions
  stamp_provenance(seu)
}

#' Diagnostic: is this clustering showing a sample/batch effect rather than
#' (or in addition to) real biology? Quantitative counterpart to
#' scrnaseq-dimensionality-reduction's visual check (R-DR-007) - rule
#' R-CLU-006.
#'
#' Neither number is a verdict on its own - a cluster dominated by one sample
#' can be a genuinely rare, sample-specific cell state, and even a modest
#' centroid separation can matter for a subtle biological signal - but both
#' are evidence worth reporting before treating a multi-sample clustering as
#' settled.
#'
#' @param seu         Clustered Seurat object.
#' @param cluster_col Metadata column with cluster assignments.
#' @param sample_col  Metadata column identifying sample/batch.
#' @param reduction   Embedding for the centroid-distance check, e.g. "umap".
#' @param dims        Which dimensions of `reduction` to use.
#' @return A list with `purity` (per-cluster data.frame) and
#'         `centroid_distances` (a `dist` object, one entry per sample pair).
check_sample_mixing <- function(seu, cluster_col = "seurat_clusters",
                                sample_col = "orig.ident",
                                reduction = "umap", dims = 1:2) {
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(sprintf(paste0(
      "Column '%s' not found in metadata - was FindClusters() run, and does ",
      "it match this column name?"), cluster_col))
  }
  samples <- seu[[sample_col, drop = TRUE]]
  if (length(unique(samples)) < 2) {
    stop("Only one sample/batch found in `", sample_col, "` - nothing to ",
         "check sample mixing against.")
  }
  if (!reduction %in% names(seu@reductions)) {
    stop(sprintf("Reduction '%s' not found - compute it first (e.g. RunUMAP()).",
                 reduction))
  }
  clusters <- seu[[cluster_col, drop = TRUE]]

  # Purity: for each cluster, the fraction of its cells coming from its
  # single most common sample. 1.0 = entirely one sample; close to that
  # sample's overall share of the dataset = well mixed.
  tab <- table(cluster = clusters, sample = samples)
  purity <- apply(tab, 1, function(x) max(x) / sum(x))
  purity_df <- data.frame(
    cluster = rownames(tab), n_cells = unname(rowSums(tab)),
    purity = unname(purity),
    dominant_sample = colnames(tab)[apply(tab, 1, which.max)],
    row.names = NULL, stringsAsFactors = FALSE)

  # Centroid distances: per-sample centroid in the embedding, pairwise
  # Euclidean distance. Samples occupying different regions of the embedding
  # (large distances) is independent evidence of the same effect.
  emb <- as.data.frame(SeuratObject::Embeddings(seu, reduction)[, dims, drop = FALSE])
  centroids <- do.call(rbind, lapply(split(emb, samples), colMeans))
  centroid_distances <- stats::dist(centroids)

  list(purity = purity_df, centroid_distances = centroid_distances)
}

#' Thin wrapper around clustree, so a model doesn't have to guess the
#' metadata column prefix (which depends on the graph's assay name).
cluster_stability_tree <- function(seu, graph_assay = SeuratObject::DefaultAssay(seu)) {
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("clustree is not installed - cannot check resolution stability.")
  }
  prefix <- paste0(graph_assay, "_snn_res.")
  cols <- grep(paste0("^", graph_assay, "_snn_res"), colnames(seu@meta.data))
  if (!length(cols)) {
    stop(sprintf("No '%s_snn_res.*' columns found - was cluster_cells() run with a resolution sweep?", graph_assay))
  }
  clustree::clustree(seu@meta.data[, cols, drop = FALSE], prefix = prefix)
}

# --- Attribution -------------------------------------------------------------

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-clustering",
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
#' FindNeighbors()/FindClusters() sequence, using the same reduction/
#' resolutions just passed to them - not only after cluster_cells().
stamp_provenance <- function(seu, reduction = NULL, resolutions = NULL) {
  prov <- list(
    pack        = PACK_CITATION$pack,
    version     = PACK_CITATION$version,
    doi         = PACK_CITATION$doi,
    reduction   = reduction %||% attr(seu, "clustering_reduction"),
    resolutions = resolutions %||% attr(seu, "clustering_resolutions"),
    run_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R           = as.character(getRversion()),
    Seurat      = as.character(utils::packageVersion("Seurat"))
  )
  existing <- attr(seu, "provenance")
  attr(seu, "provenance") <- c(existing, list(prov))
  seu
}

clustering_methods_sentence <- function(seu) {
  p <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  sprintf(paste0(
    "Cells were clustered in Seurat v%s using a shared nearest-neighbour ",
    "graph (Louvain algorithm) on the '%s' reduction, swept across ",
    "resolutions %s, following the %s workflow (%s; %s, %s, DOI: %s)."),
    p$Seurat, p$reduction, paste(range(p$resolutions), collapse = "-"),
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
