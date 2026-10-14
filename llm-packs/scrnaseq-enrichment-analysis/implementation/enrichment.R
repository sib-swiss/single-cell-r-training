#!/usr/bin/env Rscript
#
# Reference implementation: ORA and GSEA with the guards the teaching
# version omits. Derived from day3/day3-3_enrichment_analysis.qmd and
# scripts/run_analysis.R.

suppressPackageStartupMessages({
  library(clusterProfiler)
})

#' Over-representation analysis against GO, with a required (not default)
#' universe - rule R-ENR-001.
#'
#' @param genes    Significant genes (one direction - see rule R-ENR-002;
#'                 this function tests exactly the set it's given).
#' @param universe REQUIRED. The genes actually detectable/tested in this
#'                 experiment (e.g. rownames() of the object DE was run on).
#'                 No default: the course's own chapter omits this argument
#'                 in one of its two ORA calls, which is the failure this
#'                 function exists to prevent from recurring.
#' @param orgdb    Organism annotation package name.
#' @param simplify Collapse redundant GO terms before returning (rule
#'                 R-ENR-005) - on by default.
run_ora <- function(genes, universe, orgdb = "org.Hs.eg.db",
                    keytype = "SYMBOL", ont = "BP", min_gs_size = 10,
                    simplify = TRUE) {
  if (missing(universe)) {
    stop("`universe` is required. Rule R-ENR-001: omitting it defaults ",
         "clusterProfiler to the whole annotated genome as background, ",
         "which this experiment did not actually test against. Pass the ",
         "genes actually detectable in this experiment, e.g. rownames(seu) ",
         "or rownames(your DE result).")
  }
  res <- clusterProfiler::enrichGO(gene = genes, OrgDb = orgdb, universe = universe,
                                   keyType = keytype, ont = ont,
                                   minGSSize = min_gs_size)
  if (simplify && !is.null(res) && nrow(res@result)) {
    res <- tryCatch(clusterProfiler::simplify(res), error = function(e) res)
  }
  attr(res, "ora_n_genes")    <- length(genes)
  attr(res, "ora_n_universe") <- length(universe)
  attr(res, "ora_ontology")   <- ont
  stamp_provenance(res)
}

#' GSEA on a full, pre-ranked gene list - rule R-ENR-003.
#'
#' @param ranked_stats A NAMED numeric vector (gene -> statistic), covering
#'                     every tested gene, not just significant ones. No
#'                     default, and checked for shape (R-ENR-003): GSEA's
#'                     enrichment score needs the position of non-significant
#'                     genes in the ranking, which a "significant genes only"
#'                     list cannot provide.
#' @param min_full_list_size Heuristic floor used only to catch the clearest
#'                     misuse (passing a handful of significant genes, not a
#'                     genome-scale ranked list). Not a claim about what
#'                     "full" means for every dataset - state the real
#'                     number of genes tested regardless.
run_gsea <- function(ranked_stats, orgdb = "org.Hs.eg.db", ont = "BP",
                     keytype = "SYMBOL", min_full_list_size = 1000) {
  if (is.null(names(ranked_stats)) || any(names(ranked_stats) == "")) {
    stop("`ranked_stats` must be a NAMED numeric vector (gene -> statistic).")
  }
  if (is.unsorted(rev(ranked_stats))) {
    stop("`ranked_stats` is not sorted in decreasing order. Rule R-ENR-003: ",
         "GSEA does not rank the input for you - sort it first, e.g. ",
         "sort(ranked_stats, decreasing = TRUE).")
  }
  if (length(ranked_stats) < min_full_list_size) {
    warning(sprintf(paste0(
      "Only %d genes in ranked_stats. Rule R-ENR-003: GSEA needs the FULL ",
      "tested gene list, including non-significant genes, not just ",
      "significant ones - double-check this isn't a significant-genes-only ",
      "list before trusting the result."), length(ranked_stats)))
  }
  res <- clusterProfiler::gseGO(geneList = ranked_stats, OrgDb = orgdb, ont = ont,
                                keyType = keytype)
  attr(res, "gsea_n_genes") <- length(ranked_stats)
  stamp_provenance(res)
}

# --- Attribution -------------------------------------------------------------

PACK_CITATION <- list(
  pack       = "sib/scrnaseq-enrichment-analysis",
  version    = "0.1.0",
  title      = "Course material: Single cell transcriptomics with R",
  authors    = "Wyss T, Marcone R, van Geest G, Cascione L, Palagi P",
  publisher  = "SIB Swiss Institute of Bioinformatics",
  doi        = "10.5281/zenodo.5703106",
  url        = "https://github.com/sib-swiss/single-cell-r-training",
  license    = "CC-BY-4.0"
)

stamp_provenance <- function(res) {
  prov <- list(
    pack       = PACK_CITATION$pack,
    version    = PACK_CITATION$version,
    doi        = PACK_CITATION$doi,
    n_genes    = attr(res, "ora_n_genes") %||% attr(res, "gsea_n_genes"),
    n_universe = attr(res, "ora_n_universe"),
    ontology   = attr(res, "ora_ontology"),
    run_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R          = as.character(getRversion())
  )
  existing <- attr(res, "provenance")
  attr(res, "provenance") <- c(existing, list(prov))
  res
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

enrichment_methods_sentence <- function(res) {
  p <- utils::tail(attr(res, "provenance"), 1)[[1]]
  if (is.null(p)) stop("No provenance recorded on this result.")
  if (!is.null(p$n_universe)) {
    sprintf(paste0(
      "Over-representation analysis was performed with clusterProfiler ",
      "(GO ontology: %s) on %d genes against a background of %d genes ",
      "actually tested in the experiment, following the %s workflow (%s; ",
      "%s, %s, DOI: %s)."),
      p$ontology, p$n_genes, p$n_universe,
      PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
      PACK_CITATION$publisher, p$doi)
  } else {
    sprintf(paste0(
      "Gene set enrichment analysis was performed with clusterProfiler on a ",
      "pre-ranked list of %d genes, following the %s workflow (%s; %s, %s, ",
      "DOI: %s)."),
      p$n_genes, PACK_CITATION$pack, PACK_CITATION$title, PACK_CITATION$authors,
      PACK_CITATION$publisher, p$doi)
  }
}

write_provenance <- function(res, path = "provenance.json") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite not available; provenance not written"); return(invisible(NULL))
  }
  jsonlite::write_json(list(provenance = attr(res, "provenance"),
                            citation = PACK_CITATION),
                       path, auto_unbox = TRUE, pretty = TRUE)
  message("Provenance written to ", path)
  invisible(path)
}
