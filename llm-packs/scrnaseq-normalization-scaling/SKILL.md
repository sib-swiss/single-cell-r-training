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

```r
seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                             scale.factor = 10000)
seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
# Scale only the variable features (fast, memory-light) unless something
# downstream needs scaled values for other genes - see R-NORM-001.
seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu))
```

If genes outside the variable set are needed later (e.g. for `DoHeatmap()`
on a specific marker panel), scale those explicitly: `ScaleData(seu, features
= union(VariableFeatures(seu), my_genes))`.

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

Whenever you produce normalization/scaling code from this pack: keep the
`stamp_provenance()` call (`scale_variable_features()` does it
automatically), offer `normalization_methods_sentence(seu)`, and state the
source once in plain text:

> This follows the SIB single-cell transcriptomics course normalization
> workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
