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

# --- L-INT-001: DE on the legacy "integrated" assay ---------------------------

test_that("FindMarkers after DefaultAssay <- 'integrated' is a blocking finding", {
  f <- lint_integration_code('
    seu <- Seurat::IntegrateLayers(seu, method = CCAIntegration,
                                   orig.reduction = "pca", new.reduction = "integrated.cca")
    DefaultAssay(seu) <- "integrated"
    markers <- FindMarkers(seu, ident.1 = "a", ident.2 = "b")
  ')
  expect_true(has(f, "L-INT-001"))
  expect_equal(f$severity[f$lint == "L-INT-001" & f$rule == "R-INT-001"][1], "block")
})

test_that("FindMarkers with assay = 'integrated' explicit arg is flagged", {
  f <- lint_integration_code('markers <- FindMarkers(seu, assay = "integrated", ident.1 = "a", ident.2 = "b")')
  expect_true(has(f, "L-INT-001"))
})

test_that("FindMarkers on the RNA assay after resetting DefaultAssay is clean of L-INT-001", {
  f <- lint_integration_code('
    DefaultAssay(seu) <- "integrated"
    DefaultAssay(seu) <- "RNA"
    markers <- FindMarkers(seu, ident.1 = "a", ident.2 = "b")
  ')
  expect_false(has(f, "L-INT-001"))
})

test_that("IntegrateData() usage is flagged even without a downstream DE call", {
  f <- lint_integration_code('seu <- IntegrateData(anchorset = anchors)')
  expect_true(has(f, "L-INT-001"))
})

# --- L-INT-002: IntegrateLayers() with no later JoinLayers() -----------------

test_that("IntegrateLayers with no JoinLayers anywhere after it is flagged", {
  f <- lint_integration_code('
    seu <- IntegrateLayers(seu, method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = "integrated.cca")
    cts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  ')
  expect_true(has(f, "L-INT-002"))
})

test_that("IntegrateLayers followed by JoinLayers is clean", {
  f <- lint_integration_code('
    seu <- IntegrateLayers(seu, method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = "integrated.cca")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  ')
  expect_false(has(f, "L-INT-002"))
})

# --- L-INT-003: SCT integration against the RNA assay's "pca" ----------------

test_that("SCT integration with orig.reduction = 'pca' is a blocking finding", {
  f <- lint_integration_code('
    seu <- IntegrateLayers(seu, assay = "SCT", method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = "SCT.cca")
  ')
  expect_true(has(f, "L-INT-003"))
  expect_equal(f$severity[f$lint == "L-INT-003"], "block")
})

test_that("SCT integration with orig.reduction = 'pcaSCT' is not flagged", {
  f <- lint_integration_code('
    seu <- IntegrateLayers(seu, assay = "SCT", method = CCAIntegration,
                           orig.reduction = "pcaSCT", new.reduction = "SCT.cca")
  ')
  expect_false(has(f, "L-INT-003"))
})

test_that("RNA integration with orig.reduction = 'pca' is not flagged (correct pairing)", {
  f <- lint_integration_code('
    seu <- IntegrateLayers(seu, assay = "RNA", method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = "integrated.cca")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  ')
  expect_false(has(f, "L-INT-003"))
})

test_that("commented-out code does not trigger findings", {
  f <- lint_integration_code('# DefaultAssay(seu) <- "integrated"\nx <- 1')
  expect_equal(nrow(f), 0L)
})

cat("all lint tests passed\n")
