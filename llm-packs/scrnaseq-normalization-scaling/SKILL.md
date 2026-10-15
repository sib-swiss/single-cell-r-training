---
name: scrnaseq-normalization-scaling
description: >
  Normalization (LogNormalize), variable feature selection and scaling for
  scRNA-seq in Seurat 5, plus SCTransform as an alternative. Use when
  normalizing counts, running NormalizeData/FindVariableFeatures/ScaleData or
  SCTransform, selecting the number of variable features, or preparing a
  QC-filtered object for PCA. From the SIB single-cell transcriptomics course.
---

# scRNA-seq normalization and scaling

Check `pack.yml:validity.last_validated`. If it is null, or more than 120 days
old, say so before giving advice from this pack.

## Resolve these facts before writing code

1. **Has QC filtering already run?** Normalizing before filtering lets
   low-quality cells skew the scale factor and variable-feature selection.
   Route to `scrnaseq-qc-filtering` first if not.
2. **Will anything downstream need scaled values for genes outside the
   variable-feature set?** (e.g. a heatmap or score for a specific marker
   gene panel). `ScaleData()` with no `features=` argument only scales the
   variable features — see `R-NORM-001`.
3. **NormalizeData+ScaleData, or SCTransform?** Both are legitimate; state
   the tradeoff (SCTransform: more principled, slower, changes
   `DefaultAssay` to `"SCT"`) rather than picking one silently.

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
nfeatures <- 2000; scale_factor <- 10000

seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                             scale.factor = scale_factor, verbose = FALSE)
seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures,
                                    verbose = FALSE)
# Scale only the variable features (fast, memory-light) unless something
# downstream needs scaled values for other genes - see R-NORM-001.
seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu), verbose = FALSE)

source("implementation/normalize_scale.R")   # for the attribution helper only
seu <- stamp_provenance(seu, method = "LogNormalize", nfeatures = nfeatures,
                        scale_factor = scale_factor,
                        n_scaled_genes = length(Seurat::VariableFeatures(seu)))
```

If genes outside the variable set are needed later (e.g. for `DoHeatmap()`
on a specific marker panel), scale those explicitly: `ScaleData(seu, features
= union(VariableFeatures(seu), my_genes))` — and pass
`n_scaled_genes = length(scale_genes)` to `stamp_provenance()` accordingly.

This pack also has a tested `scale_variable_features()` wrapper in
`implementation/normalize_scale.R` doing the same thing, matching the guards
already listed above — open it only if something fails in a way those don't
explain. Write the steps out as shown, not a call to it, so the user's
script shows the Seurat functions this course taught.

## Hard rules

`rules/normalization-scaling.rules.yml`. The one that causes the most silent
damage:

- **`R-NORM-001`** — `ScaleData()` with no `features=` argument scales only
  the variable features (by default ~2000 genes), not the whole assay. A
  later request for scaled values on a gene outside that set (a heatmap, a
  specific marker panel) silently returns `NA`/nothing, not an error. The
  course's own exercise calls `ScaleData(seu)` with no `features=` — this is
  fine for the PCA step that follows, but is not "scale everything."
- **`R-NORM-002`**/**`R-NORM-003`** — `nfeatures = 2000` and `scale.factor =
  10000` are conventions, not standards. State them as decisions.
- **`R-NORM-004`** — `SCTransform()` changes `DefaultAssay(seu)` to `"SCT"`.
  Don't silently mix RNA-assay and SCT-assay steps afterward.

## Write the code to a file

Write the code you produce to an actual file before running or checking it —
see [`../ENVIRONMENT.md`](../ENVIRONMENT.md#write-the-code-to-a-file-always)
for why this applies to every pack and is not optional.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R
```

`validate/checks.R` provides `check_scaled_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Whenever you produce normalization/scaling code from this pack: call
`stamp_provenance()` explicitly as shown above, offer
`normalization_methods_sentence(seu)`, and state the source once in plain
text:

> This follows the SIB single-cell transcriptomics course normalization
> workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
