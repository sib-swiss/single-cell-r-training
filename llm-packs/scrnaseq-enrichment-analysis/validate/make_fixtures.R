#!/usr/bin/env Rscript
#
# Fixtures for CI. Unlike the other packs' fully synthetic gene names,
# enrichGO()/gseGO() need REAL gene symbols that resolve against real GO
# annotations - org.Hs.eg.db is a local annotation package (no network
# needed), so this stays fast and offline while still testing against real
# biology. A planted cell-cycle gene set gives the ORA test something real
# to find.

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
})

set.seed(42)

CELL_CYCLE_GENES <- c("CCNB1","CCNB2","CDK1","CCNA2","PLK1","AURKA","AURKB",
                      "BUB1","MAD2L1","CDC20","CCNE1","CCNE2","CDK2","E2F1",
                      "MCM2")

make_fixture <- function(n_universe = 800) {
  all_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  background  <- sample(setdiff(all_symbols, CELL_CYCLE_GENES), n_universe - length(CELL_CYCLE_GENES))
  universe <- union(CELL_CYCLE_GENES, background)
  list(genes = CELL_CYCLE_GENES, universe = universe)
}

#' A full ranked list for GSEA: cell-cycle genes near the top, everything
#' else at random values.
make_ranked_fixture <- function(n_total = 1500) {
  all_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  background  <- sample(setdiff(all_symbols, CELL_CYCLE_GENES), n_total - length(CELL_CYCLE_GENES))
  stats <- c(stats::runif(length(CELL_CYCLE_GENES), 3, 5),
            stats::rnorm(length(background), 0, 1))
  names(stats) <- c(CELL_CYCLE_GENES, background)
  sort(stats, decreasing = TRUE)
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fx <- make_fixture()
  saveRDS(fx, file.path(outdir, "ora_input.rds"))
  cat(sprintf("ora_input: %d genes, %d universe -> ora_input.rds\n",
              length(fx$genes), length(fx$universe)))
}
