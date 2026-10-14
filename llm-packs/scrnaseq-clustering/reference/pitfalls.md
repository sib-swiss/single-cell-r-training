# Pitfalls

---

## P-CLU-001 — A resolution presented as the correct number of clusters

**Wrong.** "Using resolution 0.3, I found the correct 9 clusters in this
dataset."

**Why it looks correct.** `FindClusters()` is deterministic for a given
resolution and random seed, and `clustree` shows the partition is "stable" at
that resolution — both make the result feel authoritative.

**What actually happened.** Source: course lecture, stated directly and
repeatedly: "clustering is subjective, there is no ground truth... you don't
know what are the true T cells and what are the true B cells in your data,
you can only have an idea." `FindClusters()` optimises a modularity function
*for the resolution it's given* — it finds a good partition at that
resolution, not the correct resolution for the dataset. A different,
equally defensible resolution gives a different, equally "optimal" (for
itself) partition.

**Detection.** Any statement that a specific resolution's output is "the
correct" or "the optimal" clustering, rather than a choice made for a stated
purpose (e.g. "broad cell types" vs. "fine subtypes").

**Right.** State the resolution as a decision tied to the granularity
wanted, show the `clustree` stability plot that informed it, and say so
explicitly — same framing as the QC pack's `R-010` for thresholds.

---

## P-CLU-002 — Clustering on the unintegrated PCA of an integrated object

**Wrong**, on a multi-sample object that has already been through
integration:

```r
seu <- Seurat::FindNeighbors(seu, dims = 1:25)   # reduction defaults to "pca"
seu <- Seurat::FindClusters(seu, resolution = 0.3)
Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.3")  # looks fine
```

**Why it looks correct.** The code runs without error, the object still has
both a `"pca"` and an `"integrated.cca"` reduction (integration adds a
reduction, it doesn't remove the original PCA), and the resulting UMAP/plot
looks like a normal clustering output — there is no obvious visual cue that
the wrong embedding was used.

**What actually happened.** `FindNeighbors()`'s `reduction=` argument
defaults to `"pca"`. On an object that has been integrated specifically to
correct for sample/batch effects, this default silently ignores that
correction and builds the neighbour graph on the pre-correction embedding —
clusters will still partly track sample identity rather than purely cell
type. The course's own exercise sets `reduction = "integrated.cca"`
explicitly, which is itself the evidence that this default is not meant to
be relied on after integration.

**Detection.** On an object carrying an integrated/harmony/cca reduction,
check that `FindNeighbors()`'s `reduction=` was set to that reduction's name,
not left to default.

**Right.**

```r
seu <- Seurat::FindNeighbors(seu, dims = 1:25, reduction = "integrated.cca")
```

---

## P-CLU-003 — Treating a coarse clustering as final with no sub-structure check

**Wrong.** Clustering once at a chosen resolution, annotating the resulting
clusters, and stopping there.

**Why it looks correct.** The clusters look clean on the UMAP and each one
has plausible marker genes — nothing in the output flags that a cluster
might be a mixture.

**What actually happened.** Source: course lecture: a cluster that looks
visually homogeneous at a coarse resolution can contain real sub-structure
that only resolves once it's analysed on its own — the recommended practice
is to subset a cluster of interest and redo scaling → PCA → UMAP →
clustering within it, to check whether it splits further.

**Right.** For any cluster the analysis depends on being homogeneous,
subset it and re-run the pipeline within it before relying on that
assumption.
