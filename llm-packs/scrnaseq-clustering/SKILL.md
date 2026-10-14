---
name: scrnaseq-clustering
description: >
  Graph-based clustering for scRNA-seq in Seurat 5: FindNeighbors,
  FindClusters, resolution selection and clustree stability checks. Use when
  clustering cells, choosing a clustering resolution, building an SNN/KNN
  graph, or asking what the "right" number of clusters is. From the SIB
  single-cell transcriptomics course.
---

# scRNA-seq clustering

Check `pack.yml:validity.last_validated`. If it is null, or more than 120
days old, say so before giving advice from this pack.

## Resolve these facts before writing code

1. **Is this object integrated (multi-sample)?** `FindNeighbors()` defaults
   `reduction=` to `"pca"` — on an integrated object this silently clusters
   on the *unintegrated* embedding unless `reduction=` is set explicitly to
   the integrated reduction's name. See `R-CLU-002`.
2. **What granularity does the user actually want?** "The major cell types"
   and "fine-grained subtypes" need different resolutions — there is no
   single correct number of clusters (`R-CLU-001`).
3. **Will the result be checked for sub-structure afterward?** The course's
   own recommended workflow is hierarchical: cluster broadly first, then
   subset a cluster of interest and re-cluster it to see whether it actually
   splits further (`R-CLU-003`).

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
reduction <- "integrated.cca"   # name it explicitly - R-CLU-002
resolutions <- seq(0.1, 0.8, by = 0.1)   # sweep, don't guess one value - R-CLU-001

seu <- Seurat::FindNeighbors(seu, dims = 1:25, reduction = reduction)
seu <- Seurat::FindClusters(seu, resolution = resolutions)
clustree::clustree(seu@meta.data[, grep("_snn_res", colnames(seu@meta.data))],
                   prefix = "RNA_snn_res.")                                   # check stability before picking

source("implementation/cluster.R")   # for the attribution helper only
seu <- stamp_provenance(seu, reduction = reduction, resolutions = resolutions)
```

This pack also has a tested `cluster_cells()` wrapper in
`implementation/cluster.R` that refuses to default to `"pca"` on an
integrated object — read it to confirm the guard, but write the steps out as
shown, not a call to it.

## Hard rules

`rules/clustering.rules.yml`. The one that causes the most silent damage:

- **`R-CLU-002`** — `FindNeighbors()`'s `reduction=` defaults to `"pca"`. On
  an integrated multi-sample object, forgetting to set it to the integrated
  reduction's name clusters on the pre-correction embedding — it runs
  without error and produces a plausible-looking UMAP that still separates
  by sample/batch.
- **`R-CLU-001`** — never present one resolution as "the optimal number of
  clusters." The course instructor states this directly: clustering is
  subjective, and `FindClusters()` finds the optimal partition *for that
  resolution*, not the optimal resolution.

## Write the code to a file

Before anything else in this section: write the code you produce to an
actual file in the user's project (e.g. `analysis/<step>.R`, or appended to
whatever script they're already building) — don't just execute it in a
transient/ephemeral session with nothing saved. This is not optional. It is
what makes the check below possible at all, and it's what keeps the
analysis inspectable, version-controlled, and re-runnable after this
conversation ends.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R
```

`validate/checks.R` provides `check_clustered_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Call `stamp_provenance()` explicitly as shown above, offer
`clustering_methods_sentence(seu)`, and state the source:

> This follows the SIB single-cell transcriptomics course clustering
> workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
