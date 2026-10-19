#!/usr/bin/env Rscript
library(testthat)

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}
source(file.path(script_dir(), "lint_code.R"))

has <- function(f, lint) lint %in% f$lint

test_that("slot= is flagged", {
  f <- lint_annotation_code('cts <- GetAssayData(seu, slot = "data")')
  expect_true(has(f, "L-ANN-001"))
})

test_that("layer= is not flagged", {
  f <- lint_annotation_code('cts <- GetAssayData(seu, layer = "data")')
  expect_false(has(f, "L-ANN-001"))
})

test_that("SingleR() with no diagnostic plot anywhere is flagged", {
  f <- lint_annotation_code('res <- SingleR::SingleR(test = m, ref = ref, labels = labs)')
  expect_true(has(f, "L-ANN-002"))
})

test_that("SingleR() WITH a diagnostic plot call is not flagged", {
  f <- lint_annotation_code('
    res <- SingleR::SingleR(test = m, ref = ref, labels = labs)
    SingleR::plotScoreHeatmap(res)
  ')
  expect_false(has(f, "L-ANN-002"))
})

test_that("manual cell-type naming with no database evidence is flagged", {
  f <- lint_annotation_code('seu <- Seurat::RenameIdents(seu, "0" = "T cell", "1" = "B cell")')
  expect_true(has(f, "L-ANN-003"))
})

test_that("manual cell-type naming WITH a marker_db lookup is not flagged", {
  f <- lint_annotation_code('
    marker_db <- readxl::read_excel(here::here("Cell_marker_Human.xlsx"))
    seu <- Seurat::RenameIdents(seu, "0" = "T cell", "1" = "B cell")
  ')
  expect_false(has(f, "L-ANN-003"))
})

test_that("manual cell-type naming WITH an explicit, code-level memory fallback is not flagged", {
  # Must be an actual code-level marker, not a comment: strip_comments()
  # removes comments before this check runs, which is intentional - the
  # rule itself requires stating the memory fallback in the reply, not only
  # in a code comment, so a comment alone should not satisfy this lint either.
  f <- lint_annotation_code('
    seu <- Seurat::RenameIdents(seu, "0" = "T cell")
    seu <- stamp_provenance(seu, marker_source = "model_memory_unverified")
  ')
  expect_false(has(f, "L-ANN-003"))
})

cat("all lint tests passed\n")
