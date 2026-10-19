#!/usr/bin/env Rscript
#
# Reference implementation: automated reference-based cell annotation.
# Derived from day3/day3-1_cell_annotation.qmd.
#
# `ref`/`labels` are passed in rather than hardcoded to a specific celldex
# dataset, so this can be tested against a small synthetic reference without
# a network call - and so a real caller can supply whichever reference suits
# their organism/tissue (rule: automated annotation needs a reference that
# actually contains the expected cell types).

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

#' Annotate cells against a reference, the SingleR way
#'
#' @param seu             Clustered Seurat object.
#' @param ref             Reference expression matrix/SummarizedExperiment
#'                        (e.g. from celldex).
#' @param labels          Reference labels, one per column of `ref`.
#' @param min_label_cells Labels assigned to fewer cells than this are
#'                        relabelled "none" (R-ANN rare-label transparency -
#'                        matches the course's own convention, informational
#'                        not silent).
annotate_with_reference <- function(seu, ref, labels, min_label_cells = 10) {
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    stop("SingleR is not installed.")
  }
  # R-ANN-001: layer=, not slot= (the default already resolves to "data").
  test_mat <- SeuratObject::GetAssayData(seu, layer = "data")

  result <- SingleR::SingleR(test = test_mat, ref = ref, labels = labels)

  tab <- table(result$labels)
  rare <- names(tab)[tab < min_label_cells]
  final_labels <- result$labels
  n_relabelled <- sum(final_labels %in% rare)
  if (length(rare)) {
    final_labels[final_labels %in% rare] <- "none"
  }

  seu$SingleR_annot <- final_labels
  attr(seu, "annotation_result")       <- result
  attr(seu, "annotation_n_labels")     <- length(unique(result$labels))
  attr(seu, "annotation_n_rare_relabelled") <- n_relabelled
  stamp_provenance(seu)
}

#' Match each cluster's top marker genes against a trusted marker-gene table
#' (e.g. downloaded from CellMarker or PanglaoDB), instead of naming a cell
#' type from recalled/memorized knowledge - rule R-ANN-007.
#'
#' Column names are explicit arguments rather than hardcoded, because a
#' downloaded database table's schema is not something to assume - verify it
#' with colnames(marker_db) first (CellMarker's own site layout changed
#' during this pack's authoring, which is exactly the kind of drift this
#' guards against).
#'
#' @param cluster_markers Data frame shaped like FindAllMarkers() output,
#'                        with at least `cluster` and `gene` columns -
#'                        already filtered to each cluster's top markers.
#' @param marker_db       Data frame with at least a cell-type column and a
#'                        gene column (e.g. CellMarker's `cell_name`/`Symbol`,
#'                        or PanglaoDB's `cell type`/`official gene symbol`).
#' @param celltype_col    Name of `marker_db`'s cell-type column.
#' @param gene_col        Name of `marker_db`'s gene column.
#' @return A data.frame, one row per cluster, with its best-matching cell
#'         type and the overlap size - a starting point to confirm against a
#'         FeaturePlot, not a label to accept unexamined.
match_clusters_to_markers <- function(cluster_markers, marker_db,
                                      celltype_col, gene_col) {
  if (!all(c("cluster", "gene") %in% names(cluster_markers))) {
    stop("`cluster_markers` needs 'cluster' and 'gene' columns (FindAllMarkers()-shaped).")
  }
  if (!all(c(celltype_col, gene_col) %in% names(marker_db))) {
    stop(sprintf(paste0(
      "`marker_db` is missing '%s' and/or '%s' - check the downloaded ",
      "table's actual column names (colnames(marker_db)) before calling ",
      "this; a database's schema is not something to assume."),
      celltype_col, gene_col))
  }
  gene_sets <- split(marker_db[[gene_col]], marker_db[[celltype_col]])
  clusters <- unique(cluster_markers$cluster)
  res <- lapply(clusters, function(cl) {
    cl_genes <- cluster_markers$gene[cluster_markers$cluster == cl]
    overlap <- vapply(gene_sets, function(g) length(intersect(cl_genes, g)), integer(1))
    best <- names(overlap)[which.max(overlap)]
    data.frame(cluster = cl, best_match = best, n_overlap = max(overlap),
              n_cluster_genes = length(cl_genes), stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

#' Both SingleR diagnostic plots in one call, so neither is forgotten
#' (rule R-ANN-002: a full annotation with no quality check is not trustworthy).
annotation_diagnostic_plots <- function(seu) {
  result <- attr(seu, "annotation_result")
  if (is.null(result)) stop("No annotation_result recorded - run annotate_with_reference() first.")
  list(
    score_heatmap     = SingleR::plotScoreHeatmap(result),
    delta_distribution = SingleR::plotDeltaDistribution(result)
  )
}

# --- Attribution -------------------------------------------------------------

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-cell-annotation",
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
#' SingleR()-based annotation, using the same label counts just computed -
#' not only after annotate_with_reference(). `marker_source` is for the
#' manual/marker-based path (rule R-ANN-007): one of "CellMarker",
#' "PanglaoDB", "connected_tool", or "model_memory_unverified" - NULL when
#' this call follows the automated SingleR path instead.
stamp_provenance <- function(seu, n_labels = NULL, n_rare_relabelled = NULL,
                             marker_source = NULL) {
  prov <- list(
    pack            = PACK_CITATION$pack,
    version         = PACK_CITATION$version,
    doi             = PACK_CITATION$doi,
    n_labels        = n_labels %||% attr(seu, "annotation_n_labels"),
    n_rare_relabelled = n_rare_relabelled %||% attr(seu, "annotation_n_rare_relabelled"),
    marker_source   = marker_source,
    run_at          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R               = as.character(getRversion()),
    Seurat          = as.character(utils::packageVersion("Seurat"))
  )
  existing <- attr(seu, "provenance")
  attr(seu, "provenance") <- c(existing, list(prov))
  seu
}

annotation_methods_sentence <- function(seu) {
  p <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  if (!is.null(p$marker_source)) {
    source_note <- if (identical(p$marker_source, "model_memory_unverified")) {
      "cell-type marker knowledge recalled from the model's training data, NOT independently verified against a database - treat these labels with caution"
    } else {
      sprintf("a reference marker table from %s", p$marker_source)
    }
    return(sprintf(paste0(
      "Clusters were manually annotated by matching marker genes against ",
      "%s, following the %s workflow (%s; %s, %s, DOI: %s)."),
      source_note, PACK_CITATION$pack, PACK_CITATION$title,
      PACK_CITATION$authors, PACK_CITATION$publisher, p$doi))
  }
  sprintf(paste0(
    "Cell types were annotated with SingleR against a reference dataset ",
    "(%d candidate labels; %d cells with a rare (<10-cell) label relabelled ",
    "'none'), following the %s workflow (%s; %s, %s, DOI: %s)."),
    p$n_labels, p$n_rare_relabelled,
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
