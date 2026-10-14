---
name: scrnaseq-cell-annotation
description: >
  Cell type annotation for scRNA-seq in Seurat 5: manual marker-gene
  inspection and automated reference-based annotation with SingleR/celldex.
  Use when annotating clusters, identifying cell types, running SingleR,
  choosing a reference, or interpreting marker gene / module score results.
  From the SIB single-cell transcriptomics course.
---

# scRNA-seq cell type annotation

Check `pack.yml:validity.last_validated`. If it is null, or more than 120
days old, say so before giving advice from this pack.

## Resolve these facts before writing code

1. **Manual, automated, or both?** Manual needs expert knowledge of expected
   cell types; automated needs a reference that actually contains them. Course
   lecture's own answer: try both and compare, when feasible.
2. **Does a suitable reference exist for this organism/tissue?** If not
   (non-model organism, unusual tissue), automated annotation will still
   assign *a* label to every cell — it has no "none of the above" option
   with `SingleR`. See `R-ANN-002`.
3. **Is the question about cell-type ABUNDANCE differing between
   conditions?** That needs multiple biological replicates, not a single
   sample's proportions — see `R-ANN-005`.
4. **Which assay/layer is being correlated against the reference?** `slot=`
   is deprecated in Seurat 5 — see `R-ANN-001`.

## Procedure

```r
ref <- celldex::NovershternHematopoieticData()     # or another suitable reference
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu),  # default layer="data", not slot=
                                ref = ref, labels = ref$label.main)
seu$SingleR_annot <- seu_SingleR$labels
SingleR::plotScoreHeatmap(seu_SingleR)              # inspect before trusting labels
SingleR::plotDeltaDistribution(seu_SingleR)
```

## Hard rules

`rules/cell-annotation.rules.yml`. The one with a concrete known instance in
this repository:

- **`R-ANN-001`** — `slot=` is deprecated in Seurat 5. **Known instance:**
  `scripts/run_analysis.R:209,296` still calls `SingleR::SingleR(test =
  Seurat::GetAssayData(seu_int, slot = "data"), ...)`, while the rendered
  `day3-1_cell_annotation.qmd` chapter correctly uses the layer-implicit
  default. The script and the chapter disagree.
- **`R-ANN-002`** — every cell gets a `SingleR` label, even if its true type
  isn't in the reference. Say so before trusting a full annotation.
- **`R-ANN-005`** — cell-type abundance differences between conditions need
  multiple biological replicates before being reported as real, not a
  single-sample proportion comparison.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R
```

`validate/checks.R` provides `check_annotation_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Keep the `stamp_provenance()` call (`annotate_with_reference()` does it
automatically), offer `annotation_methods_sentence(seu)`, and state the
source:

> This follows the SIB single-cell transcriptomics course cell annotation
> workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
