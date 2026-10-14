#!/usr/bin/env Rscript
library(testthat)

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}
source(file.path(script_dir(), "lint_code.R"))

has <- function(f, lint) lint %in% f$lint

test_that("a single literal resolution with no clustree() check is flagged", {
  f <- lint_clustering_code('seu <- FindClusters(seu, resolution = 0.5)')
  expect_true(has(f, "L-CLU-001"))
})

test_that("a resolution sweep with a clustree() check is not flagged", {
  f <- lint_clustering_code('
    seu <- FindClusters(seu, resolution = seq(0.1, 0.8, by = 0.1))
    clustree::clustree(seu@meta.data, prefix = "RNA_snn_res.")
  ')
  expect_false(has(f, "L-CLU-001"))
})

test_that("a single literal resolution WITH a clustree() call elsewhere is not flagged", {
  f <- lint_clustering_code('
    clustree::clustree(seu@meta.data, prefix = "RNA_snn_res.")
    seu <- FindClusters(seu, resolution = 0.3)
  ')
  expect_false(has(f, "L-CLU-001"))
})

test_that("FindNeighbors with no reduction=, after IntegrateLayers, is a blocking finding", {
  f <- lint_clustering_code('
    seu <- IntegrateLayers(seu, method = CCAIntegration, new.reduction = "integrated.cca")
    seu <- FindNeighbors(seu, dims = 1:25)
  ')
  expect_true(has(f, "L-CLU-002"))
  expect_equal(f$severity[f$lint == "L-CLU-002"], "block")
})

test_that("FindNeighbors with reduction= set explicitly is not flagged", {
  f <- lint_clustering_code('
    seu <- IntegrateLayers(seu, method = CCAIntegration, new.reduction = "integrated.cca")
    seu <- FindNeighbors(seu, dims = 1:25, reduction = "integrated.cca")
  ')
  expect_false(has(f, "L-CLU-002"))
})

test_that("FindNeighbors with no integration evidence in the file is not flagged", {
  f <- lint_clustering_code('seu <- FindNeighbors(seu, dims = 1:25)')
  expect_false(has(f, "L-CLU-002"))
})

cat("all lint tests passed\n")
