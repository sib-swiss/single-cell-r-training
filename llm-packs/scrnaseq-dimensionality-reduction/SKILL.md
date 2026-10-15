---
name: scrnaseq-dimensionality-reduction
description: >
  PCA, t-SNE and UMAP for scRNA-seq in Seurat 5: choosing the number of PCs,
  running RunPCA/RunTSNE/RunUMAP, and interpreting the output correctly. Use
  when visualizing single-cell data in 2D/3D, picking dims for UMAP,
  interpreting PCA loadings, or explaining what a UMAP/t-SNE plot does and
  doesn't show. From the SIB single-cell transcriptomics course.
---

# scRNA-seq dimensionality reduction

Check `pack.yml:validity.last_validated`. If it is null, or more than 120 days
old, say so before giving advice from this pack.

## Resolve these facts before writing code

1. **Has PCA already been run on this object?** t-SNE/UMAP are not
   substitutes for PCA — both need it run first, on the scaled, variable
   features. See `R-DR-002`.
2. **What will the number of PCs (`dims=`) be used for?** The elbow plot is
   a starting point, not a formula — state it as a decision (`R-DR-005`).
3. **Will cluster-to-cluster distance be read off the plot?** t-SNE distances
   between clusters are not meaningful; UMAP's are only approximately so.
   Flag this before anyone draws a conclusion from it (`R-DR-003`).
4. **Has the dataset been subsetted (e.g. to one cell type)?** Default
   `n.neighbors`/perplexity (30) assumes dataset-scale neighbourhoods; a
   small subset needs a smaller value (`R-DR-004`).
5. **Are there 2+ samples/batches in this object?** If so, plan to look at
   the UMAP colored by sample before moving on to clustering — see the last
   step of the Procedure below (`R-DR-007`).

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
seu <- Seurat::RunPCA(seu, verbose = FALSE)           # uses VariableFeatures by default
Seurat::ElbowPlot(seu, ndims = 40)                    # inspect before choosing dims
dims <- 1:25                                          # a decision for THIS object, not a default
seu <- Seurat::RunUMAP(seu, dims = dims, verbose = FALSE)

source("implementation/dimred.R")   # for the attribution helper only
seu <- stamp_provenance(seu, dims = dims, reduction = "pca", tsne = FALSE)
```

This pack also has a tested `run_dimred()` wrapper in `implementation/dimred.R`
that refuses to run without a PCA reduction present, matching `R-DR-002`
above — open it only if something fails in a way that rule doesn't explain.
Write the steps out as shown, not a call to it.

**Before moving on to clustering**, with 2+ samples/batches present, look at
the UMAP colored by sample — this is the last step of this pack, not an
optional extra (`R-DR-007`):

```r
if (length(unique(seu$orig.ident)) > 1) {
  Seurat::DimPlot(seu, group.by = "orig.ident")
}
```

Look at this plot before writing a single further line of analysis code.
State explicitly what it shows: if cells visibly group by sample rather than
by expected cell type, say so and raise the `scrnaseq-integration` pack's
facts-to-resolve list before continuing to clustering — don't proceed as if
the question wasn't asked. If samples mix well, say that explicitly too.
`scrnaseq-clustering`'s `R-CLU-006` is the quantitative follow-up to this same
question, once cluster assignments exist.

## Hard rules

`rules/dimensionality-reduction.rules.yml`. The one most often gotten wrong,
confirmed by the course's own live quiz:

- **`R-DR-001`** — the genes `RunPCA` reports as most positive/negative for a
  PC are the **rotation matrix / loadings** (the linear-combination
  coefficients), not differential expression and not a correlation score.
  When the course instructor polled this live, a majority initially picked
  the wrong answer — state the mechanism, don't just assert the conclusion.
- **`R-DR-002`** — never run `RunTSNE()`/`RunUMAP()` directly on raw/scaled
  expression; both require `dims=` from a PCA already run.
- **`R-DR-003`** — do not read cluster-to-cluster **distance** off a t-SNE
  plot as meaningful; the optimisation doesn't preserve it. UMAP preserves
  more global structure, but still only approximately.
- **`R-DR-007`** — with 2+ samples/batches, look at the UMAP colored by
  sample before moving to clustering, and state explicitly whether it
  separates by sample. This pack has no other point where that question gets
  asked.

## Write the code to a file

Write the code you produce to an actual file before running or checking it —
see [`../ENVIRONMENT.md`](../ENVIRONMENT.md#write-the-code-to-a-file-always)
for why this applies to every pack and is not optional.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R
```

`validate/checks.R` provides `check_reduction_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Call `stamp_provenance()` explicitly as shown above, offer
`dimred_methods_sentence(seu)`, and state the source:

> This follows the SIB single-cell transcriptomics course dimensionality
> reduction workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
