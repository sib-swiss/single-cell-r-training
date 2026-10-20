# Pitfalls

---

## P-NORM-001 — Assuming `ScaleData()` scales every gene

**Wrong.**

```r
seu <- Seurat::ScaleData(seu)
# ... much later ...
Seurat::DoHeatmap(seu, features = my_favourite_20_marker_genes)
# some rows of the heatmap are blank/NA
```

**Why it looks correct.** `ScaleData(seu)` runs without error, `seu` gets a
populated `scale.data` layer, and nothing in the console output suggests
anything was skipped.

**What actually happened.** `ScaleData()`'s default `features` argument is
`VariableFeatures(object = seu)`, not every gene in the assay. Genes outside
that ~2000-gene variable set simply have no entry in `scale.data`. This is
invisible until something asks for one of those genes specifically — a
heatmap of a curated marker panel, a module score computed from scaled
values — at which point it returns empty/`NA` rather than an error that
explains why.

**Detection.**

```r
stopifnot(all(my_genes %in% rownames(SeuratObject::LayerData(seu, layer = "scale.data"))))
```

**Right.** For the standard PCA workflow this is exactly correct — PCA only
ever uses variable features, so scaling only them is the efficient, intended
behaviour, and the course's own exercise call (`ScaleData(seu)`, no
`features=`) is fine for that purpose. It stops being fine the moment
anything downstream needs scaled values for a specific gene outside that
set: pass `features` explicitly in that case.

---

## P-NORM-002 — `scale.factor`/`nfeatures` presented as universal standards

**Wrong.** "I normalized with the standard scale factor of 10,000 and
selected the standard 2,000 variable features."

**Why it looks correct.** These are literally the values in Seurat's own
documentation and in this course's exercise, used by default in nearly every
published analysis.

**What actually happened.** `scale.factor` only needs to be large relative to
typical per-cell library size so that normalized values aren't awkwardly
small; 10,000 is a convention inherited from bulk RNA-seq's CPM, not a
property of the data being analysed. `nfeatures` trades off capturing more
biological variation (higher) against computational cost and noise from
rarely-informative genes (lower) — the right value depends on how
heterogeneous the sample is, not on a fixed convention.

**Right.** State both as decisions with a stated rationale, the same
pattern as the QC pack's thresholds (`R-010`).

---

## P-NORM-003 — Mixing RNA-assay and SCT-assay steps without noticing

**Wrong.**

```r
seu <- Seurat::SCTransform(seu)
# ... several steps later, in a different function/script ...
seu <- Seurat::ScaleData(seu)   # intending to scale the RNA assay
```

**Why it looks correct.** Every call runs without error. `ScaleData()`
doesn't know or care that the user meant the RNA assay.

**What actually happened.** `SCTransform()` sets `DefaultAssay(seu)` to
`"SCT"`. Every subsequent assay-implicit call (`ScaleData()`, `RunPCA()`,
etc. with no explicit `assay=`) now operates on the SCT assay, silently,
until something downstream that specifically expects the RNA assay (e.g.
this course's QC pack's guards, or a DE pipeline built around raw/LogNormalize
counts) gets data it wasn't written for.

**Detection.** After `SCTransform()`, check `DefaultAssay(seu)` before
trusting what any subsequent assay-implicit call operated on.

**Right.** Reset explicitly: `DefaultAssay(seu) <- "RNA"` immediately after
`SCTransform()` if anything downstream expects the RNA assay, and say so.
