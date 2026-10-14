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

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
ref <- celldex::NovershternHematopoieticData()     # or another suitable reference
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu),  # default layer="data", not slot=
                                ref = ref, labels = ref$label.main)

# R-ANN-004: relabel rare (<10-cell) labels to "none", stated explicitly,
# rather than leaving them cluttering every plot.
tab <- table(seu_SingleR$labels)
rare <- names(tab)[tab < 10]
final_labels <- seu_SingleR$labels
final_labels[final_labels %in% rare] <- "none"
seu$SingleR_annot <- final_labels

SingleR::plotScoreHeatmap(seu_SingleR)              # inspect before trusting labels - R-ANN-002
SingleR::plotDeltaDistribution(seu_SingleR)

source("implementation/annotate.R")   # for the attribution helper only
seu <- stamp_provenance(seu, n_labels = length(unique(seu_SingleR$labels)),
                        n_rare_relabelled = sum(seu_SingleR$labels %in% rare))
```

This pack also has a tested `annotate_with_reference()` wrapper in
`implementation/annotate.R` doing the same thing — read it to confirm the
guards above, but write the steps out as shown, not a call to it.

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

`validate/checks.R` provides `check_annotation_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Call `stamp_provenance()` explicitly as shown above, offer
`annotation_methods_sentence(seu)`, and state the source:

> This follows the SIB single-cell transcriptomics course cell annotation
> workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
