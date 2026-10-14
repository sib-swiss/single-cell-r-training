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

# --- L-DE-001: cross-condition per-cell test with no pseudobulk step --------

test_that("FindMarkers across a condition with no pseudobulk step is blocking", {
  f <- lint_de_code('
    seu <- SetIdent(seu, value = "type")
    de <- FindMarkers(seu, ident.1 = "disease", ident.2 = "healthy", group.by = "type")
  ')
  expect_true(has(f, "L-DE-001"))
  expect_equal(f$severity[f$lint == "L-DE-001"], "block")
})

test_that("FindMarkers across a condition WITH pseudobulk elsewhere in the file is clean of L-DE-001", {
  f <- lint_de_code('
    bulk <- AggregateExpression(seu, group.by = "orig.ident", return.seurat = TRUE)
    de <- FindMarkers(seu, ident.1 = "disease", ident.2 = "healthy", group.by = "type")
  ')
  expect_false(has(f, "L-DE-001"))
})

test_that("FindAllMarkers for cluster discovery (no group.by condition) is not flagged", {
  f <- lint_de_code('markers <- FindAllMarkers(seu, min.pct = 0.25, only.pos = TRUE)')
  expect_false(has(f, "L-DE-001"))
})

# --- L-DE-002: pseudobulk pooling across replicates within a condition ------

test_that("AggregateExpression grouped only by condition is flagged", {
  f <- lint_de_code('bulk <- AggregateExpression(seu, group.by = "type", return.seurat = TRUE)')
  expect_true(has(f, "L-DE-002"))
})

test_that("AggregateExpression grouped by sample is not flagged", {
  f <- lint_de_code('bulk <- AggregateExpression(seu, group.by = "orig.ident", return.seurat = TRUE)')
  expect_false(has(f, "L-DE-002"))
})

# --- L-DE-003: significant gene list with no diagnostics --------------------

test_that("filtering to p_val_adj < 0.05 with no diagnostic plot is flagged", {
  f <- lint_de_code('sig <- subset(de, de$p_val_adj < 0.05)')
  expect_true(has(f, "L-DE-003"))
})

test_that("filtering to significance WITH de_diagnostic_plots() call is clean", {
  f <- lint_de_code('
    sig <- subset(de, de$p_val_adj < 0.05)
    plots <- de_diagnostic_plots(de)
  ')
  expect_false(has(f, "L-DE-003"))
})

# --- L-DE-005: repeated pairwise FindMarkers with no omnibus alternative ----

test_that("two separate FindMarkers calls with no omnibus alternative are flagged", {
  f <- lint_de_code('
    de_young_old <- FindMarkers(seu, ident.1 = "young", ident.2 = "old")
    de_young_mid <- FindMarkers(seu, ident.1 = "young", ident.2 = "middle")
  ')
  expect_true(has(f, "L-DE-005"))
})

test_that("two FindMarkers calls alongside anova_de() are not flagged", {
  f <- lint_de_code('
    omnibus <- anova_de(seu, sample_col = "orig.ident", group_col = "type")
    de_young_old <- FindMarkers(seu, ident.1 = "young", ident.2 = "old")
    de_young_mid <- FindMarkers(seu, ident.1 = "young", ident.2 = "middle")
  ')
  expect_false(has(f, "L-DE-005"))
})

test_that("a single FindMarkers call is not flagged by L-DE-005", {
  f <- lint_de_code('de <- FindMarkers(seu, ident.1 = "a", ident.2 = "b")')
  expect_false(has(f, "L-DE-005"))
})

test_that("commented-out code does not trigger findings", {
  f <- lint_de_code('# de <- FindMarkers(seu, ident.1 = "a", ident.2 = "b", group.by = "type")\nx <- 1')
  expect_equal(nrow(f), 0L)
})

cat("all lint tests passed\n")
