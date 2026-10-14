# Pitfalls

---

## P-ANN-001 — `slot=` instead of `layer=` feeding SingleR

**Wrong.**

```r
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu, slot = "data"),
                                ref = ref, labels = ref$label.main)
```

**Why it looks correct.** It runs — `slot=` is deprecated, not removed, in
current Seurat, and this is a *real, live instance* in this repository:
`scripts/run_analysis.R:209` and `:296` both use exactly this pattern.

**What actually happened.** Same failure family as the QC pack's `R-004`:
`slot=` is the Seurat 4 argument name; `layer=` is current. The rendered
`day3-1_cell_annotation.qmd` chapter already uses the layer-implicit default
(`Seurat::GetAssayData(seu)`, no `slot=`/`layer=` at all) — the script and
the chapter disagree, and the script is what a learner is more likely to
copy wholesale.

**Right.** `Seurat::GetAssayData(seu)` (defaults to `layer = "data"`), or
explicit `layer = "data"`.

---

## P-ANN-002 — Trusting a complete `SingleR` annotation with no quality check

**Wrong.** Running `SingleR()`, assigning `seu$SingleR_annot <-
seu_SingleR$labels`, and treating every cell's label as settled.

**Why it looks correct.** `SingleR()` returns a label for every single cell
— there is no missing data, no obvious gap to prompt suspicion.

**What actually happened.** Source: course lecture, stated directly:
`SingleR` has no "rejection option." If the sample contains a cell type that
simply isn't in the reference panel, those cells still get assigned *some*
label — the best-correlating one available, which is wrong. Nothing in the
basic output distinguishes a confident, correct label from a forced,
wrong one.

**Detection.** Annotation reported with no accompanying
`plotScoreHeatmap()`/`plotDeltaDistribution()` inspection.

**Right.** Generate both diagnostic plots, look for cells/clusters with weak
or ambiguous score separation between labels, and say explicitly that
`SingleR` will not flag a genuinely novel cell type on its own — that has to
come from the diagnostics, or from checking the reference's coverage before
trusting the annotation.

---

## P-ANN-003 — Single marker gene mRNA level treated as a surface-identity readout

**Wrong.** "This cell doesn't express `CD4`, so it isn't CD4-positive."

**Why it looks correct.** `CD4` is the textbook marker for the cell surface
protein of the same name — using the gene's mRNA level as a proxy for the
protein's presence is standard practice.

**What actually happened.** Source: course lecture, with a concrete example
from a CITE-seq dataset (surface protein measured directly alongside mRNA):
cells with clearly present CD4 surface protein nonetheless showed low CD4
mRNA. The mRNA-to-surface-protein relationship is not always tight, and
single-cell RNA-seq only measures mRNA.

**Right.** Use several marker genes for a cell type call, not one, and treat
a single gene's low expression as inconclusive rather than exclusionary.

---

## P-ANN-004 — Cell-type abundance differences reported from one sample per condition

**Wrong.** "Condition B has almost no cluster 3 cells compared to condition
A — this cell type is depleted in condition B," from one sample of each.

**Why it looks correct.** The proportions in the data genuinely differ, and
sometimes strikingly so — the course's own example shows a population nearly
disappearing between conditions.

**What actually happened.** Source: course lecture, a direct caution: "it's
not always the case that the frequencies that you have in your dataset...
is really a good image of what you had in your actual sample... you always
have to validate if your differences are actually true." A single sample's
cell-type proportions are noisy (capture efficiency, dissociation bias, and
sampling all affect them) and are not automatically a reliable estimate of
the true biological proportion. The course's own example that *did* trust
such a difference only did so because it was consistent across 3 biological
replicates per condition, and was independently confirmed by flow cytometry.

**Detection.** A cell-type abundance comparison between conditions with no
stated replicate count, or n=1 per condition.

**Right.** State the number of biological replicates behind the comparison,
and flag that an orthogonal validation (e.g. flow cytometry) strengthens the
claim considerably — don't present a single-sample proportion difference as
settled biology.
