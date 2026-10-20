#!/usr/bin/env Rscript
library(testthat)

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}
source(file.path(script_dir(), "lint_code.R"))

has <- function(f, lint) lint %in% f$lint

test_that("enrichGO() with no universe= is a blocking finding", {
  f <- lint_enrichment_code('
    res <- clusterProfiler::enrichGO(genes, "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  ')
  expect_true(has(f, "L-ENR-001"))
  expect_equal(f$severity[f$lint == "L-ENR-001"], "block")
})

test_that("enrichGO() with universe= is not flagged", {
  f <- lint_enrichment_code('
    res <- clusterProfiler::enrichGO(genes, "org.Hs.eg.db",
                                     universe = rownames(proB),
                                     keyType = "SYMBOL", ont = "BP")
  ')
  expect_false(has(f, "L-ENR-001"))
})

test_that("GSEA() with no sort/order call earlier is a blocking finding", {
  f <- lint_enrichment_code('res <- clusterProfiler::GSEA(geneList = stats)')
  expect_true(has(f, "L-ENR-002"))
})

test_that("GSEA() preceded by a sort() call is not flagged", {
  f <- lint_enrichment_code('
    stats <- sort(stats, decreasing = TRUE)
    res <- clusterProfiler::GSEA(geneList = stats)
  ')
  expect_false(has(f, "L-ENR-002"))
})

test_that("enrichGO() result accessed with no simplify() is flagged", {
  f <- lint_enrichment_code('
    res <- clusterProfiler::enrichGO(genes, "org.Hs.eg.db", universe = u)
    View(res@result)
  ')
  expect_true(has(f, "L-ENR-003"))
})

test_that("enrichGO() result accessed WITH simplify() is not flagged", {
  f <- lint_enrichment_code('
    res <- clusterProfiler::enrichGO(genes, "org.Hs.eg.db", universe = u)
    res <- clusterProfiler::simplify(res)
    View(res@result)
  ')
  expect_false(has(f, "L-ENR-003"))
})

cat("all lint tests passed\n")
