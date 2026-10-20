#!/usr/bin/env Rscript
# Tests for the static linter. Runs without Seurat.
#   Rscript validate/test-lint.R

library(testthat)

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}
source(file.path(script_dir(), "lint_code.R"))

has <- function(f, lint) lint %in% f$lint

# --- the known-bad case from the course material ------------------------------

bad_hardcoded_mt <- '
seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
'

test_that("hardcoded ^MT- with no species check is a blocking finding", {
  f <- lint_qc_code(bad_hardcoded_mt)
  expect_true(has(f, "L-001"))
  expect_equal(f$severity[f$lint == "L-001"], "block")
})

test_that("unguarded PercentageFeatureSet is a blocking finding", {
  f <- lint_qc_code(bad_hardcoded_mt)
  expect_true(has(f, "L-002"))
})

test_that("threshold without a plot is flagged", {
  f <- lint_qc_code(bad_hardcoded_mt)
  expect_true(has(f, "L-005"))
})

# --- the corrected version should be clean ------------------------------------

good <- '
# species resolved by the caller
mito_pattern <- switch(species, human = "^MT-", mouse = "^mt-")
mito_genes <- grep(mito_pattern, rownames(seu), value = TRUE)
stopifnot("no mitochondrial genes matched" = length(mito_genes) > 0)
seu <- Seurat::PercentageFeatureSet(seu, pattern = mito_pattern,
                                    col.name = "percent.mito")
Seurat::VlnPlot(seu, features = c("nFeature_RNA", "percent.mito"),
                group.by = "orig.ident")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
'

test_that("the guarded, species-aware version is clean", {
  f <- lint_qc_code(good)
  expect_equal(nrow(f), 0L)
})

# --- targeted checks ----------------------------------------------------------

test_that("slot= is flagged as a Seurat 5 problem", {
  f <- lint_qc_code('cts <- GetAssayData(seu, assay = "RNA", slot = "counts")')
  expect_true(has(f, "L-003"))
})

test_that("layer= is not flagged", {
  f <- lint_qc_code('cts <- GetAssayData(seu, assay = "RNA", layer = "counts")')
  expect_false(has(f, "L-003"))
})

test_that("over-broad haemoglobin pattern is flagged", {
  f <- lint_qc_code('seu <- PercentageFeatureSet(seu, pattern = "^HB", col.name = "pg")')
  expect_true(has(f, "L-004"))
})

test_that("the course haemoglobin pattern is not flagged as over-broad", {
  f <- lint_qc_code('seu <- PercentageFeatureSet(seu, pattern = "^HB[^(P)]", col.name = "pg")')
  expect_false(has(f, "L-004"))
})

test_that("commented-out code does not trigger findings", {
  f <- lint_qc_code('# seu <- PercentageFeatureSet(seu, pattern = "^MT-")\nx <- 1')
  expect_equal(nrow(f), 0L)
})

test_that("a mouse-correct pattern is not flagged by L-001", {
  f <- lint_qc_code('
    mito_genes <- grep("^mt-", rownames(seu), value = TRUE)
    stopifnot(length(mito_genes) > 0)
    seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mito")
    VlnPlot(seu, "percent.mito")
  ')
  expect_equal(nrow(f), 0L)
})

# --- regression: whole-file species scan produced a false negative -------------
# scripts/run_analysis.R mentions "Human" ~180 lines below the mito pattern, in
# an unrelated celldex call. That must not satisfy the species check.

test_that("a distant, unrelated species mention does not suppress L-001", {
  far <- c('seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")',
           rep("x <- 1", 40),
           'ref <- celldex::HumanPrimaryCellAtlasData()')
  f <- lint_qc_code(far)
  expect_true(has(f, "L-001"))
})

test_that("a nearby species determination does suppress L-001", {
  near <- c('species <- "human"',
            'mito_genes <- grep("^MT-", rownames(seu), value = TRUE)',
            'stopifnot(length(mito_genes) > 0)',
            'seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")',
            'VlnPlot(seu, "percent.mito")')
  f <- lint_qc_code(near)
  expect_false(has(f, "L-001"))
})


# --- L-005 fires on literals, not on parameterised thresholds -----------------
# Regression: implementation/qc_filter.R:filter_cells() takes thresholds as
# arguments and leaves plotting to its caller. That is the shape R-010 wants.

test_that("a parameterised threshold is not flagged by L-005", {
  f <- lint_qc_code('
    filter_cells <- function(seu, min_features, max_mito) {
      seu[, seu$nFeature_RNA > min_features & seu$percent.mito < max_mito]
    }')
  expect_false(has(f, "L-005"))
})

test_that("a hardcoded literal threshold with no plot is flagged by L-005", {
  f <- lint_qc_code('seu <- subset(seu, subset = percent.mito < 8)')
  expect_true(has(f, "L-005"))
})

# --- L-006 (R-004): subset() after split() with no JoinLayers() --------------

test_that("subset() after split() with no JoinLayers() is flagged", {
  f <- lint_qc_code('
    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
    seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
  ')
  expect_true(has(f, "L-006"))
})

test_that("subset() after split() AND JoinLayers() is not flagged", {
  f <- lint_qc_code('
    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
    seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
  ')
  expect_false(has(f, "L-006"))
})

test_that("subset() with no split() anywhere is not flagged by L-006", {
  f <- lint_qc_code('seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)')
  expect_false(has(f, "L-006"))
})

# --- L-007 (R-006): MALAT1 removed with no stated reason ----------------------

test_that("MALAT1 removed with no nearby justification is flagged", {
  f <- lint_qc_code('seu <- seu[rownames(seu) != "MALAT1", ]')
  expect_true(has(f, "L-007"))
})

test_that("MALAT1 removed with a stated QC-evidence reason nearby is not flagged", {
  f <- lint_qc_code('
    # percent.mito is also elevated in these cells, consistent with stress
    seu <- seu[rownames(seu) != "MALAT1", ]
  ')
  expect_false(has(f, "L-007"))
})

test_that("MALAT1 merely mentioned, not removed, is not flagged", {
  f <- lint_qc_code('cat("Most expressed gene is MALAT1, consistent with poly-A selection.\\n")')
  expect_false(has(f, "L-007"))
})

# --- L-008 (R-007): raw (not filtered) CellRanger matrix ----------------------

test_that("reading raw_feature_bc_matrix is flagged", {
  f <- lint_qc_code('sparse_matrix <- Seurat::Read10X(data.dir = "outs/raw_feature_bc_matrix")')
  expect_true(has(f, "L-008"))
})

test_that("reading filtered_feature_bc_matrix is not flagged", {
  f <- lint_qc_code('sparse_matrix <- Seurat::Read10X(data.dir = "outs/filtered_feature_bc_matrix")')
  expect_false(has(f, "L-008"))
})

# --- L-009 (R-009): PercentageFeatureSet after NormalizeData/SCTransform -----

test_that("PercentageFeatureSet after NormalizeData is a blocking finding", {
  f <- lint_qc_code('
    seu <- Seurat::NormalizeData(seu)
    seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")
  ')
  expect_true(has(f, "L-009"))
  expect_equal(f$severity[f$lint == "L-009"], "block")
})

test_that("PercentageFeatureSet before NormalizeData is not flagged", {
  f <- lint_qc_code('
    seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")
    seu <- Seurat::NormalizeData(seu)
  ')
  expect_false(has(f, "L-009"))
})

# --- L-010 (R-014): percent.ribo used as a filtering threshold ----------------

test_that("percent.ribo used in a filter threshold is flagged", {
  f <- lint_qc_code('seu <- subset(seu, subset = percent.ribo > 5)')
  expect_true(has(f, "L-010"))
})

test_that("percent.ribo merely reported (VlnPlot) is not flagged", {
  f <- lint_qc_code('Seurat::VlnPlot(seu, features = "percent.ribo")')
  expect_false(has(f, "L-010"))
})

cat("all lint tests passed\n")
