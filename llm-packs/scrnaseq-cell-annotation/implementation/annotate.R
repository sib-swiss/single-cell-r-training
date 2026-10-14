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

stamp_provenance <- function(seu) {
  prov <- list(
    pack            = PACK_CITATION$pack,
    version         = PACK_CITATION$version,
    doi             = PACK_CITATION$doi,
    n_labels        = attr(seu, "annotation_n_labels"),
    n_rare_relabelled = attr(seu, "annotation_n_rare_relabelled"),
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
