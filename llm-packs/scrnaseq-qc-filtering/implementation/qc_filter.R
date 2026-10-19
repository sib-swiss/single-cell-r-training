#!/usr/bin/env Rscript
#
# Reference implementation: scRNA-seq QC and cell filtering, 10x + Seurat 5.
#
# This is the artifact an LLM should ADAPT rather than reimplement. It is
# exercised by validate/run_all.R in CI against synthetic human, mouse and
# Drosophila fixtures, so its behaviour on the failure modes in
# rules/qc-filtering.rules.yml is a tested property, not a claim.
#
# Derived from day1/day1-3_analysis_tools_qc.qmd, with the guards that the
# teaching version omits for readability.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

# --- Species-specific feature patterns ---------------------------------------
#
# The single most important table in this pack. Drosophila patterns are taken
# from scripts/drosophila_rp_mt.R in this repo - evidence that the course
# already knows these differ per species, but had not encoded it anywhere a
# learner copying scripts/run_analysis.R would encounter.

QC_PATTERNS <- list(
  human = list(
    mito   = "^MT-",
    ribo   = "^RP[SL]",
    globin = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBZ")
  ),
  mouse = list(
    mito   = "^mt-",
    ribo   = "^Rp[sl]",
    globin = c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Hbb-y","Hbb-bh1")
  ),
  drosophila = list(
    mito   = "^mt:",
    ribo   = "^Rp[SL]",
    globin = character(0)          # no haemoglobin orthologue
  )
)

#' Add per-cell QC metrics, failing loudly on an empty pattern match
#'
#' @param seu     A Seurat object.
#' @param species One of names(QC_PATTERNS). No default - see rule R-001.
#' @return The object with percent.mito, percent.ribo, percent.globin added.
add_qc_metrics <- function(seu, species) {
  if (missing(species)) {
    stop("`species` is required. Guessing it is rule R-001: the wrong ",
         "mitochondrial pattern matches zero genes, silently yields ",
         "percent.mito == 0 for every cell, and disables filtering.")
  }
  species <- match.arg(species, names(QC_PATTERNS))
  p <- QC_PATTERNS[[species]]

  # R-009: QC metrics must be computed on raw counts, before normalisation.
  # PercentageFeatureSet() reads the assay's current "data" layer, which
  # CreateSeuratObject() initialises identical to "counts" - so non-integer
  # values there mean NormalizeData()/SCTransform() has already run. Sampled
  # rather than checked in full, since this only needs to catch "the whole
  # layer was rescaled", not find a single stray value.
  assay_name <- SeuratObject::DefaultAssay(seu)
  layers <- SeuratObject::Layers(seu[[assay_name]])
  if ("data" %in% layers) {
    probe <- as.numeric(SeuratObject::LayerData(seu, assay = assay_name, layer = "data")[
      seq_len(min(50, nrow(seu))), seq_len(min(20, ncol(seu))), drop = FALSE])
    if (length(probe) && any(probe != 0) && !isTRUE(all.equal(probe, round(probe)))) {
      stop("The '", assay_name, "' assay's 'data' layer does not look like raw ",
           "counts (non-integer values found). Rule R-009: QC metrics must be ",
           "computed before normalisation. Call add_qc_metrics() before running ",
           "NormalizeData or SCTransform, or recompute metrics explicitly from ",
           "the 'counts' layer.")
    }
  }

  # R-002: assert the match is non-empty BEFORE trusting the result.
  mito_genes <- grep(p$mito, rownames(seu), value = TRUE)
  if (!length(mito_genes)) {
    stop(sprintf(
      paste0("No features matched the %s mitochondrial pattern '%s'.\n",
             "  Likely causes: wrong species; features are Ensembl IDs rather ",
             "than symbols;\n  a reference that prefixes gene names; or the ",
             "default assay is not RNA.\n  First 5 features present: %s"),
      species, p$mito, paste(utils::head(rownames(seu), 5), collapse = ", ")))
  }
  seu <- Seurat::PercentageFeatureSet(seu, pattern = p$mito, col.name = "percent.mito")

  ribo_genes <- grep(p$ribo, rownames(seu), value = TRUE)
  if (!length(ribo_genes)) {
    warning(sprintf("No features matched the ribosomal pattern '%s'.", p$ribo))
    seu$percent.ribo <- NA_real_
  } else {
    seu <- Seurat::PercentageFeatureSet(seu, pattern = p$ribo, col.name = "percent.ribo")
  }

  # R-003: explicit gene list, not a regex. "^HB" would pull in HBEGF/HBS1L.
  hb <- intersect(p$globin, rownames(seu))
  if (!length(hb)) {
    seu$percent.globin <- NA_real_
  } else {
    seu <- Seurat::PercentageFeatureSet(seu, features = hb, col.name = "percent.globin")
  }

  attr(seu, "qc_species") <- species
  seu
}

#' Summarise QC distributions per sample
#'
#' R-005: thresholds chosen from a pooled distribution hide between-sample
#' differences in capture quality.
qc_summary <- function(seu, group.by = "orig.ident") {
  md <- seu@meta.data
  grp <- if (group.by %in% colnames(md)) md[[group.by]] else rep("all", nrow(md))
  metrics <- intersect(c("nCount_RNA","nFeature_RNA","percent.mito",
                         "percent.ribo","percent.globin"), colnames(md))
  do.call(rbind, lapply(split(seq_len(nrow(md)), grp), function(idx) {
    row <- data.frame(sample = as.character(grp[idx][1]), n_cells = length(idx))
    for (m in metrics) {
      q <- stats::quantile(md[[m]][idx], c(0.05, 0.5, 0.95), na.rm = TRUE)
      row[[paste0(m, "_p05")]] <- unname(q[1])
      row[[paste0(m, "_med")]] <- unname(q[2])
      row[[paste0(m, "_p95")]] <- unname(q[3])
    }
    row
  }))
}

#' Apply QC filters, reporting the effect
#'
#' Thresholds are REQUIRED arguments with no defaults. That is deliberate and is
#' rule R-010: a default here would be copied as a standard, and there is no
#' such thing. The values used in the course - min_features 200, max_features
#' 5000, max_mito 8 - come from Caron et al. 2020 for paediatric bone marrow.
filter_cells <- function(seu, min_features, max_features, max_mito,
                         max_removed_frac = 0.40, verbose = TRUE) {
  stopifnot(!missing(min_features), !missing(max_features), !missing(max_mito))

  # R-013: filtering after integration invalidates the learned correction.
  if (length(seu@reductions) &&
      any(grepl("integrated|harmony|cca|mnn", names(seu@reductions), ignore.case = TRUE))) {
    stop("This object already carries an integrated reduction (",
         paste(names(seu@reductions), collapse = ", "),
         "). QC filtering belongs before integration - rule R-013.")
  }
  if (is.null(seu$percent.mito)) {
    stop("percent.mito is absent; call add_qc_metrics() first.")
  }
  # R-004/P-012: subset()-style filtering on an unjoined v5 object does not
  # error - it can silently operate within a single layer rather than across
  # every sample.
  assay_name <- SeuratObject::DefaultAssay(seu)
  layers <- tryCatch(SeuratObject::Layers(seu[[assay_name]]), error = function(e) character(0))
  n_count_layers <- sum(grepl("^counts", layers))
  if (n_count_layers > 1) {
    stop(sprintf(paste0(
      "%d unjoined 'counts' layers found on assay '%s' (%s). Call ",
      "SeuratObject::JoinLayers() before filtering - rule R-004/P-012."),
      n_count_layers, assay_name, paste(grep("^counts", layers, value = TRUE), collapse = ", ")))
  }
  # R-001 runtime post-condition.
  if (all(seu$percent.mito == 0, na.rm = TRUE)) {
    stop("percent.mito is 0 for every cell. This is nearly always a pattern/",
         "species mismatch rather than unusually clean data - rule R-001.")
  }

  n_before <- table(seu$orig.ident)
  keep <- seu$nFeature_RNA > min_features &
          seu$nFeature_RNA < max_features &
          seu$percent.mito < max_mito
  seu_f <- seu[, keep]
  n_after <- table(seu_f$orig.ident)

  report <- data.frame(
    sample  = names(n_before),
    before  = as.integer(n_before),
    after   = as.integer(n_after[names(n_before)]),
    stringsAsFactors = FALSE)
  report$after[is.na(report$after)] <- 0L
  report$removed_frac <- 1 - report$after / report$before

  overall <- 1 - ncol(seu_f) / ncol(seu)
  if (verbose) {
    message(sprintf("Filtering: %d -> %d cells (%.1f%% removed)",
                    ncol(seu), ncol(seu_f), 100 * overall))
    print(report, row.names = FALSE)
  }
  # R-011
  if (overall > max_removed_frac) {
    stop(sprintf(paste0("Filtering removed %.1f%% of cells, above the %.0f%% ",
                        "guard. Re-inspect the thresholds against the ",
                        "distributions before proceeding - rule R-011."),
                 100 * overall, 100 * max_removed_frac))
  }
  attr(seu_f, "qc_report") <- report
  attr(seu_f, "qc_thresholds") <- list(min_features = min_features,
                                       max_features = max_features,
                                       max_mito = max_mito)
  stamp_provenance(seu_f)
}

# --- Attribution -------------------------------------------------------------
#
# The attribution problem for an LLM pack is that provenance is stripped at
# every hop: pack -> model -> generated code -> the user's script -> the
# methods section. By the time a paper is written, the researcher may not know
# the pack existed - they asked an assistant a question and got an answer.
#
# So attribution has to ride with the ARTIFACT rather than the human's memory,
# and has to be something the researcher actively wants. Hence: a machine
# -readable stamp on the object, and a ready-to-paste methods sentence.

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-qc-filtering",
  version    = "0.1.0",
  title      = "Course material: Single cell transcriptomics with R",
  authors    = "Wyss T, Marcone R, van Geest G, Cascione L, Palagi P",
  publisher  = "SIB Swiss Institute of Bioinformatics",
  doi        = "10.5281/zenodo.5703106",
  url        = "https://github.com/sib-swiss/single-cell-r-training",
  license    = "CC-BY-4.0"
)

#' Record which pack, version and parameters produced this object
#'
#' Stored as an attribute and, where possible, in the Seurat object's own
#' command log, so it survives saveRDS() and travels into whatever the
#' researcher deposits alongside their paper.
#' @param species    Species used for the QC patterns. Optional - falls back
#'                   to attr(seu, "qc_species") if add_qc_metrics() set it.
#'                   Pass it explicitly when calling this after writing the
#'                   PercentageFeatureSet()/subset() steps inline rather than
#'                   through add_qc_metrics()/filter_cells().
#' @param thresholds Named list(min_features=, max_features=, max_mito=).
#'                   Same fallback as species.
stamp_provenance <- function(seu, species = NULL, thresholds = NULL) {
  prov <- list(
    pack        = PACK_CITATION$pack,
    version     = PACK_CITATION$version,
    doi         = PACK_CITATION$doi,
    species     = species %||% attr(seu, "qc_species"),
    thresholds  = thresholds %||% attr(seu, "qc_thresholds"),
    n_cells     = ncol(seu),
    run_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R           = as.character(getRversion()),
    Seurat      = as.character(utils::packageVersion("Seurat"))
  )
  existing <- attr(seu, "provenance")
  attr(seu, "provenance") <- c(existing, list(prov))
  seu
}

#' A methods sentence the researcher can paste into a manuscript
#'
#' Deliberately useful rather than merely dutiful. A pre-written, accurate,
#' parameter-filled sentence is something people want to copy; an obligation
#' to cite is something they forget.
qc_methods_sentence <- function(seu) {
  p <- utils::tail(attr(seu, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this object.")
  t <- p$thresholds
  sprintf(paste0(
    "Quality control and cell filtering were performed in Seurat v%s ",
    "following the %s workflow (%s; %s, %s, DOI: %s). ",
    "Cells were retained with %d < nFeature_RNA < %d and mitochondrial ",
    "transcript content < %g%% (%s mitochondrial gene set), leaving %d cells."),
    p$Seurat, PACK_CITATION$pack, PACK_CITATION$title,
    PACK_CITATION$authors, PACK_CITATION$publisher, p$doi,
    t$min_features, t$max_features, t$max_mito,
    p$species %||% "species-specific", p$n_cells)
}

# base R gained %||% in 4.4.0; define it only if absent
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Write a provenance sidecar next to the analysis outputs
#'
#' JSON so that it is greppable and machine-readable - this is what makes
#' downstream detection of pack use possible at all (see llm-packs/README.md).
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
