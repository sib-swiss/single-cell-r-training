#!/usr/bin/env Rscript
library(testthat)

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}
source(file.path(script_dir(), "lint_code.R"))

has <- function(f, lint) lint %in% f$lint

test_that("bare ScaleData() with a DoHeatmap() call is flagged", {
  f <- lint_norm_code('
    seu <- ScaleData(seu)
    DoHeatmap(seu, features = marker_panel)
  ')
  expect_true(has(f, "L-NORM-001"))
})

test_that("ScaleData() with explicit features= is not flagged", {
  f <- lint_norm_code('
    seu <- ScaleData(seu, features = union(VariableFeatures(seu), marker_panel))
    DoHeatmap(seu, features = marker_panel)
  ')
  expect_false(has(f, "L-NORM-001"))
})

test_that("bare ScaleData() with no DoHeatmap() in the file is not flagged", {
  f <- lint_norm_code('seu <- ScaleData(seu)')
  expect_false(has(f, "L-NORM-001"))
})

test_that("SCTransform with no DefaultAssay reset is flagged", {
  f <- lint_norm_code('seu <- SCTransform(seu)')
  expect_true(has(f, "L-NORM-002"))
})

test_that("SCTransform followed by a DefaultAssay reset is not flagged", {
  f <- lint_norm_code('
    seu <- SCTransform(seu)
    DefaultAssay(seu) <- "RNA"
  ')
  expect_false(has(f, "L-NORM-002"))
})

cat("all lint tests passed\n")
