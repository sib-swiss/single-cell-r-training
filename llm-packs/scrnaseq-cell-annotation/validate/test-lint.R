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

cat("all lint tests passed\n")
