#!/usr/bin/env Rscript
library(testthat)

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}
source(file.path(script_dir(), "lint_code.R"))

has <- function(f, lint) lint %in% f$lint

test_that("RunUMAP with no PCA evidence anywhere is a blocking finding", {
  f <- lint_dimred_code('seu <- RunUMAP(seu, dims = 1:25)')
  expect_true(has(f, "L-DR-001"))
  expect_equal(f$severity[f$lint == "L-DR-001"], "block")
})

test_that("RunUMAP after an actual RunPCA() call is not flagged", {
  f <- lint_dimred_code('
    seu <- RunPCA(seu)
    seu <- RunUMAP(seu, dims = 1:25)
  ')
  expect_false(has(f, "L-DR-001"))
})

test_that("RunUMAP after a runtime reductions-existence check is not flagged", {
  f <- lint_dimred_code('
    if (!"pca" %in% names(seu@reductions)) {
      stop("no pca reduction found")
    }
    seu <- RunUMAP(seu, dims = 1:25)
  ')
  expect_false(has(f, "L-DR-001"))
})

test_that("RunUMAP naming a non-pca reduction is not flagged", {
  f <- lint_dimred_code('seu <- RunUMAP(seu, dims = 1:25, reduction = "integrated.cca")')
  expect_false(has(f, "L-DR-001"))
})

test_that("a dims literal with no ElbowPlot anywhere is flagged", {
  f <- lint_dimred_code('seu <- RunUMAP(seu, dims = 1:25)')
  expect_true(has(f, "L-DR-002"))
})

test_that("a dims literal WITH an ElbowPlot call is not flagged", {
  f <- lint_dimred_code('
    ElbowPlot(seu, ndims = 40)
    seu <- RunUMAP(seu, dims = 1:25)
  ')
  expect_false(has(f, "L-DR-002"))
})

cat("all lint tests passed\n")
