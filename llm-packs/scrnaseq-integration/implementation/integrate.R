#!/usr/bin/env Rscript
#
# Reference implementation: multi-sample integration, 10x + Seurat 5.
#
# This is the artifact an LLM should ADAPT rather than reimplement. It is
# exercised by validate/run_all.R in CI against synthetic fixtures, so its
# behaviour on the failure modes in rules/integration.rules.yml is a tested
# property, not a claim.
#
# Derived from day2/day2-3_integration.qmd, with the guards that the teaching
# version omits for readability, plus a diagnostic for R-INT-003 the chapter
# doesn't carry a function for.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

#' Integrate multiple samples/batches with layer-based CCA
#'
#' @param seu          A Seurat object with 2+ samples and QC already applied.
#' @param assay        "RNA" or "SCT". No default - see rule R-INT-005.
#' @param batch_var    metadata column identifying samples/batches.
#' @param orig_reduction PCA to use as the integration starting point. No
#'                     default for the SCT assay (R-INT-005): the single
#'                     most damaging silent failure in this pack is reusing
#'                     the RNA assay's PCA here.
#' @param new_reduction Name for the resulting integrated reduction.
#' @param method       Integration method; defaults to Seurat::CCAIntegration,
#'                     the only one taught in this course.
#' @return The object with layers re-joined and the integrated reduction added.
integrate_samples <- function(seu, assay, batch_var = "orig.ident",
                              orig_reduction = NULL, new_reduction = NULL,
                              method = Seurat::CCAIntegration) {
  if (missing(assay)) {
    stop("`assay` is required ('RNA' or 'SCT'). Guessing it risks rule ",
         "R-INT-005: integrating the SCT assay against an RNA-assay PCA ",
         "integrates on the wrong basis without erroring.")
  }
  assay <- match.arg(assay, c("RNA", "SCT"))

  # R-INT-002 (partial, automatable half): nothing to integrate.
  batches <- seu[[batch_var, drop = TRUE]]
  if (length(unique(batches)) < 2) {
    stop("Only one batch/sample found in `", batch_var, "`. There is nothing ",
         "to integrate - proceed directly to dimensionality reduction.")
  }

  if (is.null(orig_reduction)) {
    if (assay == "SCT") {
      stop("`orig_reduction` is required when assay = 'SCT' (rule R-INT-005). ",
           "Pass the name of a PCA computed ON THE SCT ASSAY (e.g. run ",
           "RunPCA(seu, reduction.name = 'pcaSCT') first and pass ",
           "orig_reduction = 'pcaSCT'). Reusing the RNA assay's 'pca' here ",
           "integrates on the wrong basis without erroring.")
    }
    orig_reduction <- "pca"
  }
  if (!orig_reduction %in% names(seu@reductions)) {
    stop(sprintf(
      "Reduction '%s' not found on the object. Run RunPCA on the '%s' assay ",
      "before integrating it.", orig_reduction, assay))
  }
  # R-INT-005 runtime post-condition: Seurat's DimReduc objects record which
  # assay they were computed on (`assay.used`). Checking it directly is a
  # stronger guard than naming convention alone - it catches a mismatch even
  # when a reduction has a misleading name.
  used <- tryCatch(Seurat::DefaultAssay(seu@reductions[[orig_reduction]]),
                   error = function(e) NA_character_)
  if (!is.na(used) && !identical(used, assay)) {
    stop(sprintf(paste0(
      "Reduction '%s' was computed on assay '%s', but you are integrating ",
      "assay '%s'. Rule R-INT-005: a PCA computed on the wrong assay ",
      "produces anchors in a space that does not reflect the data being ",
      "corrected. Run RunPCA on the '%s' assay and pass that reduction ",
      "instead."), orig_reduction, used, assay, assay))
  }
  if (is.null(new_reduction)) {
    new_reduction <- if (assay == "SCT") "SCT.cca" else "integrated.cca"
  }

  # The course's own exercise converts "SCT" to a plain Assay5 under its own
  # name before splitting it - splitting the legacy SCTAssay class directly
  # produces a "subscript out of bounds" error from a scale.data layer-naming
  # mismatch. The RNA assay is already a plain Assay5 in Seurat 5 and needs no
  # such conversion.
  integrate_assay <- assay
  if (assay == "SCT") {
    seu[["SCT_v5"]] <- as(object = seu[["SCT"]], Class = "Assay5")
    integrate_assay <- "SCT_v5"
  }

  seu[[integrate_assay]] <- split(seu[[integrate_assay]], f = batches)

  seu <- Seurat::IntegrateLayers(
    object = seu, assay = integrate_assay, method = method,
    orig.reduction = orig_reduction, new.reduction = new_reduction,
    features = if (assay == "SCT") Seurat::VariableFeatures(seu, assay = "SCT") else NULL,
    verbose = FALSE)

  # R-INT-004: re-join before anything downstream reads a single matrix.
  seu[[integrate_assay]] <- SeuratObject::JoinLayers(seu[[integrate_assay]])

  attr(seu, "integration_method") <- "CCA"
  attr(seu, "integration_assay")  <- assay
  attr(seu, "integration_batch_var") <- batch_var
  attr(seu, "integration_reduction") <- new_reduction
  stamp_provenance(seu)
}

#' Diagnostic for R-INT-003: does a biological grouping still separate
#' post-integration, or has it been integrated away?
#'
#' Not proof either way - a non-significant result is also consistent with a
#' small effect size or few cells - but a clear positive result on a known
#' covariate is reassuring, and a clear absence where one was expected is
#' worth stopping to investigate before trusting the integrated embedding for
#' anything beyond visualisation.
#'
#' @param seu       Integrated Seurat object.
#' @param reduction The integrated reduction, e.g. "integrated.cca".
#' @param group_by  Metadata column carrying the biological comparison of
#'                  interest (NOT the batch variable that was corrected for).
#' @param dims      Which dimensions of the reduction to test.
check_integration_preserves_signal <- function(seu, reduction, group_by,
                                               dims = 1:5) {
  if (!reduction %in% names(seu@reductions)) {
    stop(sprintf("Reduction '%s' not found.", reduction))
  }
  emb <- SeuratObject::Embeddings(seu, reduction)[, dims, drop = FALSE]
  grp <- seu[[group_by, drop = TRUE]]
  if (length(unique(grp)) < 2) {
    stop("`group_by` has fewer than 2 levels - nothing to test separation for.")
  }
  pvals <- apply(emb, 2, function(x) {
    tryCatch(stats::kruskal.test(x, grp)$p.value, error = function(e) NA_real_)
  })
  res <- data.frame(dimension = colnames(emb), p_value = unname(pvals),
                    separates = unname(pvals) < 0.05,
                    stringsAsFactors = FALSE)
  if (!any(res$separates, na.rm = TRUE)) {
    message("No tested dimension of '", reduction, "' separates by '", group_by,
            "'. If a real biological difference was expected here, this is ",
            "consistent with over-integration (rule R-INT-003) - investigate ",
            "before treating the integrated embedding as trustworthy for ",
            "anything beyond batch-effect removal.")
  }
  res
}

# --- Attribution -------------------------------------------------------------
# Same mechanism as scrnaseq-qc-filtering/implementation/qc_filter.R - see that
# file's comment for the rationale.

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-integration",
  version    = "0.1.0",
  title      = "Course material: Single cell transcriptomics with R",
  authors    = "Wyss T, Marcone R, van Geest G, Cascione L, Palagi P",
  publisher  = "SIB Swiss Institute of Bioinformatics",
  doi        = "10.5281/zenodo.5703106",
  url        = "https://github.com/sib-swiss/single-cell-r-training",
  license    = "CC-BY-4.0"
)

#' Explicit arguments let this be called directly after an inline
#' split()/IntegrateLayers()/JoinLayers() sequence, using the same
#' assay/batch_var/reduction just used - not only after integrate_samples().
stamp_provenance <- function(seu, method = NULL, assay = NULL,
                             batch_var = NULL, reduction = NULL) {
  prov <- list(
    pack       = PACK_CITATION$pack,
    version    = PACK_CITATION$version,
    doi        = PACK_CITATION$doi,
    method     = method %||% attr(seu, "integration_method") %||% "CCA",
    assay      = assay %||% attr(seu, "integration_assay"),
    batch_var  = batch_var %||% attr(seu, "integration_batch_var"),
    reduction  = reduction %||% attr(seu, "integration_reduction"),
    n_cells    = ncol(seu),
    run_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R          = as.character(getRversion()),
    Seurat     = as.character(utils::packageVersion("Seurat"))
  )
  existing <- attr(seu, "provenance")
  attr(seu, "provenance") <- c(existing, list(prov))
  seu
}

#' A methods sentence the researcher can paste into a manuscript
integration_methods_sentence <- function(seu) {
  p <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  sprintf(paste0(
    "Samples were integrated in Seurat v%s using layer-based %s integration ",
    "(assay: %s, batch variable: %s, reduction: %s) following the %s ",
    "workflow (%s; %s, %s, DOI: %s), on %d cells."),
    p$Seurat, p$method, p$assay, p$batch_var, p$reduction,
    PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
    PACK_CITATION$publisher, p$doi, p$n_cells)
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

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
