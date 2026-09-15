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
  seu_f
}
