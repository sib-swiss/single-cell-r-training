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

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
seu <- Seurat::RunPCA(seu, verbose = FALSE)           # uses VariableFeatures by default
Seurat::ElbowPlot(seu, ndims = 40)                    # inspect before choosing dims
seu <- Seurat::RunUMAP(seu, dims = 1:25)               # dims is a decision, not a default
```

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

`validate/checks.R` provides `check_reduction_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Keep the `stamp_provenance()` call (`run_umap_pipeline()` does it
automatically), offer `dimred_methods_sentence(seu)`, and state the source:

> This follows the SIB single-cell transcriptomics course dimensionality
> reduction workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
