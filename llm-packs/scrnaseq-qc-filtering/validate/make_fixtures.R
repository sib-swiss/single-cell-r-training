#!/usr/bin/env Rscript
#
# Synthetic fixtures for CI. Deliberately NOT derived from course_data:
#   - no download, so CI stays fast and works on a fork
#   - deterministic, so expected outputs are exact
#   - tiny (300 cells x ~200 genes)
#   - the mouse and Drosophila fixtures exist to make the R-001 silent-failure
#     mode reproducible, which real human data cannot do
#
# Realism is only needed in the gene NAMES and the rough shape of the count
# distribution - not in the biology.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

MITO <- list(
  human      = c("MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP6","MT-CO3",
                 "MT-ND4","MT-ND5","MT-CYB"),
  mouse      = c("mt-Nd1","mt-Nd2","mt-Co1","mt-Co2","mt-Atp6","mt-Co3",
                 "mt-Nd4","mt-Nd5","mt-Cytb"),
  drosophila = c("mt:ND1","mt:ND2","mt:CoI","mt:CoII","mt:ATPase6","mt:CoIII",
                 "mt:Cyt-b")
)
RIBO <- list(
  human      = c("RPS3","RPS6","RPS19","RPL5","RPL13","RPL23"),
  mouse      = c("Rps3","Rps6","Rps19","Rpl5","Rpl13","Rpl23"),
  drosophila = c("RpS3","RpS6","RpL5","RpL13")
)
GLOBIN <- list(
  human      = c("HBA1","HBA2","HBB","HBG1"),
  mouse      = c("Hba-a1","Hbb-bs","Hbb-bt"),
  drosophila = character(0)
)
# Decoys for rule R-003: match "^HB" but are not haemoglobin.
HB_DECOY <- list(human = c("HBEGF","HBS1L","HBP1"), mouse = character(0),
                 drosophila = character(0))

make_fixture <- function(species, n_cells = 300, n_other = 180, n_samples = 3) {
  mito <- MITO[[species]]; ribo <- RIBO[[species]]
  glob <- GLOBIN[[species]]; decoy <- HB_DECOY[[species]]
  other <- sprintf("GENE%03d", seq_len(n_other))
  genes <- c(mito, ribo, glob, decoy, other)

  # Negative-binomial-ish counts; mito genes get a per-cell multiplier so that
  # percent.mito spans a realistic range with a tail of high-mito cells.
  mito_load <- c(rgamma(round(n_cells * 0.85), shape = 2, rate = 2.5),
                 rgamma(n_cells - round(n_cells * 0.85), shape = 9, rate = 1.2))
  mito_load <- sample(mito_load)

  m <- matrix(0L, nrow = length(genes), ncol = n_cells,
              dimnames = list(genes, sprintf("CELL%04d", seq_len(n_cells))))
  depth <- round(rlnorm(n_cells, meanlog = log(3000), sdlog = 0.5))
  for (j in seq_len(n_cells)) {
    w <- rep(1, length(genes))
    w[genes %in% mito] <- 3 * mito_load[j]
    w[genes %in% ribo] <- 6
    if (length(glob)) w[genes %in% glob] <- if (j %% 25 == 0) 40 else 0.2
    m[, j] <- rmultinom(1, size = depth[j], prob = w / sum(w))[, 1]
  }

  # A few deliberately poor cells so filtering has something to remove.
  bad <- sample(n_cells, round(n_cells * 0.06))
  m[, bad] <- round(m[, bad] * 0.05)

  seu <- Seurat::CreateSeuratObject(
    counts = as(m, "CsparseMatrix"), project = species,
    min.cells = 0, min.features = 0)
  seu$orig.ident <- factor(sprintf("%s-S%d", species,
                                   rep_len(seq_len(n_samples), n_cells)))
  Seurat::Idents(seu) <- seu$orig.ident
  seu
}

if (!interactive()) {
  outdir <- file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fixtures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (sp in names(MITO)) {
    seu <- make_fixture(sp)
    saveRDS(seu, file.path(outdir, sprintf("%s.rds", sp)))
    cat(sprintf("%-11s %d genes x %d cells -> %s.rds\n",
                sp, nrow(seu), ncol(seu), sp))
  }
}
