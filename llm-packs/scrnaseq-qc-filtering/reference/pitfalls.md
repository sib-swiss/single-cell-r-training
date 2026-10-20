# Pitfalls

Each entry: the wrong output, **why it looks correct**, what actually happened,
how to detect it, the right version.

Written as worked failures rather than advice because models pattern-match to
labelled anti-examples far more reliably than to prose warnings. "Don't rely on
default thresholds" is unfalsifiable and changes nothing; "a <8% mito cutoff on
cardiac tissue deletes cardiomyocytes" changes behaviour.

---

## P-001 — Mitochondrial pattern does not match the species

**Wrong**, on a mouse dataset:

```r
seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
```

**Why it looks correct.** The code is idiomatic and is what every tutorial —
including `scripts/run_analysis.R` in this repository — shows. It runs without
warning. `VlnPlot(seu, "percent.mito")` produces a plot. `subset()` succeeds.
Nothing anywhere says "0 features matched".

**What actually happened.** Mouse mitochondrial genes are `mt-Nd1`, `mt-Co1`,
`mt-Atp6`. `"^MT-"` matched nothing. `PercentageFeatureSet` returned 0 for every
cell. The filter retained 100% of cells. The object now carries a `percent.mito`
column that looks like QC was performed, and it will be described in a methods
section.

**Species patterns** (`implementation/qc_filter.R:QC_PATTERNS`):

| Species | Mitochondrial | Ribosomal |
|---|---|---|
| Human | `^MT-` | `^RP[SL]` |
| Mouse | `^mt-` | `^Rp[sl]` |
| Drosophila | `^mt:` | `^Rp[SL]` |

The Drosophila patterns come from `scripts/drosophila_rp_mt.R` in this
repository. The course already knew these differ per species — it just wasn't
encoded anywhere a learner copying `scripts/run_analysis.R` would meet it.

**Detection.**

```r
stopifnot("no mitochondrial genes matched" = length(grep(pattern, rownames(seu))) > 0)
# post hoc:
max(seu$percent.mito, na.rm = TRUE) > 0
```

**Same failure, other causes:** features are Ensembl IDs (`ENSG...`); a custom
reference prefixes every gene (`GRCh38-MT-ND1`); the default assay is `ADT` in a
CITE-seq object; a `SoupX`-corrected matrix reimported without gene symbols.

---

## P-002 — A published threshold copied as a default

**Wrong.**

> "I've applied the standard QC thresholds: `nFeature_RNA` between 200 and 5000,
> `percent.mito` below 8%."

**Why it looks correct.** These are the values in `day1-3_analysis_tools_qc.qmd`,
taken from Caron et al. 2020. They are well justified *for paediatric bone
marrow mononuclear cells on 10x 3' v2*.

**What actually happened.** There is no standard. On cardiac tissue an 8%
mitochondrial ceiling discards cardiomyocytes, which are genuinely
mitochondria-rich — the most abundant cell type in the sample is deleted and
what remains looks clean. On a shallow run the `>200` feature floor removes most
of the dataset. On brown adipose, kidney proximal tubule or skeletal muscle the
same problem recurs. Tumour tissue is a further course-lecture-sourced example:
mitochondrial gene expression in tumour cells is commonly higher than in most
other tissues, so a threshold tuned for bone marrow can be too strict there too.

**Detection.** Plot first, always. If the cutoff sits on the shoulder of a mode
rather than in a valley, it is cutting through a cell population rather than
between populations.

**Right.** Show `VlnPlot` and `FeatureScatter` split by sample, state that the
value is a decision, name its provenance, and say what would change it.

---

## P-003 — Erythroid cells discarded as contamination

**Wrong.** Filtering out high-`percent.globin` cells as "red blood cell
contamination".

**Why it looks correct.** Haemoglobin transcripts are a recognised contaminant
in many tissues, and these cells form a visually distinct outlier population.

**What actually happened.** In bone marrow — which is what the course dataset
is — those are erythroid progenitors, an expected lineage. In PBMMC-2 they
appear as high globin / low ribosomal, negatively correlated, which is the
signature of erythroid differentiation, not lysis. Removing them removes a
result.

**Right.** Establish the tissue before treating globin as contamination. Carry
them through dimensionality reduction and see where they land.

---

## P-004 — Seurat 5 layer semantics

**Wrong**, on a multi-sample v5 object:

```r
cts <- GetAssayData(seu, assay = "RNA", slot = "counts")
```

**What actually happened.** Two distinct problems. `slot=` is deprecated in
favour of `layer=`. And a v5 object built from several samples keeps counts in
per-sample layers (`counts.PBMMC-1`, `counts.PBMMC-2`, …) until `JoinLayers()`
is called — so this returns one sample's matrix, or errors, depending on
version. Any per-cell statistic computed from it silently covers a subset.

**Live instance in this repository.** `scripts/run_analysis.R:47` still uses
`slot = "counts"` inside `most_expressed_boxplot()`, while the corresponding
chapter `day1-3_analysis_tools_qc.qmd` has been updated to `layer =`. The script
and the rendered material disagree, and the script is what learners copy.

**Right.**

```r
seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
cts <- SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts")
stopifnot(ncol(cts) == ncol(seu))
```

---

## P-005 — Count-based filtering presented as doublet removal

**Wrong.** `nFeature_RNA < 5000` described as removing doublets.

**What actually happened.** It removes *some* heterotypic doublets and *some*
large healthy cells, and misses homotypic doublets entirely — two B cells look
like one B cell with somewhat more RNA. Any reported doublet rate is meaningless.

**Right.** Run a dedicated caller per sample, before merging and before
integration. `DoubletFinder` 2.0.6 is in this repo's `renv.lock`;
`scDblFinder` is generally preferred but is **not** currently in the lockfile
(see `pack.yml:depends_on`). Say explicitly that the count filter is a crude
complement, not a substitute.

---

## P-006 — `min.cells` applied before merging samples

**Wrong.**

```r
seu <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 100)
```

run separately per sample, then merged.

**What actually happened.** The retained gene set becomes sample-dependent: a
gene expressed in 3 cells of sample A but 2 of sample B is kept in one and
dropped in the other. After merging, that gene reads as zero in sample B, which
is a technical artefact indistinguishable from biology in downstream DE.

**Right.** Apply `min.cells` once, after merging — or not at all, and filter
genes explicitly with a documented rule.

---

## P-007 — Ambient RNA not corrected

**Wrong.** Treating a QC-filtered object as free of background contamination.

**Why it looks correct.** The object has been filtered for low-quality cells
and doublets. It looks clean: no implausible cells remain, every remaining
barcode passed the thresholds.

**What actually happened.** Every droplet, including real cells, contains some
ambient (free-floating, lysed-cell) RNA mixed into its count. This "soup" is
highest for whatever is most highly expressed in the sample overall, and it is
not addressed by filtering cells on `nFeature_RNA`/`percent.mito` — those
metrics operate per cell, not per gene, and a contaminated cell still looks
like a normal, healthy cell on those metrics. A gene that is truly restricted
to one cell type can appear at a low, "leaky" level everywhere.

**Detection.** Look for genes that are highly expressed in one expected
cell-type cluster and present at a low but non-zero level in clusters where
they have no plausible biological role (e.g. haemoglobin genes appearing at a
low level in T cells). This can only really be assessed after clustering, not
during per-cell filtering.

**Right.** State explicitly whether ambient correction has been applied.
`SoupX` is the standard tool for this, but is **not** currently in this
course's `renv.lock` (see `pack.yml:depends_on`) — if a user asks for ambient
RNA correction, say that it requires a package not in the course container
rather than silently skipping the step or substituting a cell-level filter
that does not address it.

---

## P-008 — QC run after integration

**Wrong.** Loading an already-integrated object and applying `filter_cells()`
or a manual `subset()` to remove low-quality cells "now that samples are
combined."

**Why it looks correct.** The combined object is convenient to work with —
one object, one set of reductions — and filtering it feels like a natural
next cleanup step.

**What actually happened.** `IntegrateLayers()` learns a correction from the
cells present at the time it ran. Removing cells afterward does not undo or
update that correction; the reduction it produced was computed against a
cell set that no longer exists in the filtered object. This is rule `R-013`,
already enforced at runtime by `filter_cells()` (it checks
`seu@reductions` for an integrated/harmony/cca reduction and refuses to run),
but is worth stating as a worked failure because it is an easy trap: the
natural, linear-feeling order is "combine samples, then clean up," when QC
and doublet removal must happen per sample, before merging and before
integration.

**Detection.**

```r
stopifnot(!any(grepl("integrated|harmony|cca",
                     names(seu@reductions), ignore.case = TRUE)))
```

**Right.** Run QC and doublet removal on each sample separately (or on the
merged-but-not-integrated object), before `IntegrateLayers()` is ever called.

---

## P-009 — A global threshold across heterogeneous samples

**Wrong.** Computing one `nFeature_RNA`/`percent.mito` threshold from the
pooled distribution across all samples in a multi-sample object, without
checking whether samples differ.

**Why it looks correct.** `VlnPlot(seu, "percent.mito")` with no `group.by`
produces one tidy violin, and one threshold is simpler to state than several.

**What actually happened.** If one sample had worse dissociation, a longer
processing delay, or simply a different capture efficiency, its distribution
sits somewhere else on the same axis. A single pooled threshold then removes
cells preferentially from the worst sample — composition differences between
samples become partly technical (how many cells survived filtering) rather
than purely biological, and that technical difference is easy to mistake for
a real effect in whatever comparison the samples were collected for.

**Detection.** Plot QC metrics split by `orig.ident` (`VlnPlot(...,
group.by = "orig.ident")`) before choosing a threshold, not after. If the
per-sample distributions visibly differ, that is itself worth stating,
whether or not the final decision is to use one threshold or several.

**Right.** State explicitly whether thresholds were set globally or per
sample, and show the split-by-sample plot that informed the decision — this
is rule `R-005`.

---

## P-010 — MALAT1 (or another highly-expressed gene) removed by habit

**Wrong.**

```r
seu <- seu[rownames(seu) != "MALAT1", ]
```

presented as routine cleanup, with no reference to anything else in the QC
output.

**Why it looks correct.** MALAT1 is very highly expressed in essentially
every 10x poly-A-enriched dataset, a lot of publicly visible analysis code
removes it unconditionally, and "remove the one gene that dominates every
plot" sounds like reasonable hygiene.

**What actually happened.** MALAT1 can carry real biological relevance (see
`reference/citations.md#c-04`); removing it is not free. The course's own QC
chapter computes `most_expressed_boxplot()`, sees MALAT1 at the top, checks
the rest of the QC picture for any other sign of stressed or dying cells,
finds none, and explicitly **leaves it in** — the decision follows from the
evidence, it is not a default step applied regardless of what the rest of QC
shows.

**Detection.** Any removal of a specific named gene with no stated connection
to the rest of the QC metrics for that dataset.

**Right.** Report the most highly expressed genes, connect the decision to
remove (or keep) any of them to the rest of the QC evidence (e.g. "percent.mito
is also elevated in the same cells, consistent with stress" vs. "no other QC
metric is unusual for these cells, so it is left in"), and say which was the
case here.

---

## P-011 — Raw vs filtered CellRanger matrix confusion

**Wrong.**

```r
sparse_matrix <- Seurat::Read10X(data.dir = "outs/raw_feature_bc_matrix")
seu <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells = 3, min.features = 100)
# ... proceed with this pack's thresholds unchanged
```

**Why it looks correct.** `Read10X()` and `CreateSeuratObject()` both run
without error on either directory — the raw matrix is a perfectly valid
sparse matrix, just a much bigger one (hundreds of thousands of barcodes
instead of a few thousand cells).

**What actually happened.** `cellranger count` already performs its own
cell-calling step on the way to producing `filtered_feature_bc_matrix/`: it
orders barcodes by UMI count and looks for the characteristic sharp drop
("knee point") that separates real cells from empty droplets, and
`filtered_feature_bc_matrix/` is the result of that step. `raw_feature_bc_matrix/`
is every barcode cellranger saw, overwhelmingly empty droplets and ambient
RNA, with a tiny fraction of real cells buried in it.
This course's exercises and this pack's implementation assume the
`filtered_feature_bc_matrix/` (post-knee-point) input — applying
`min.features = 100` or this pack's `nFeature_RNA`/`percent.mito` thresholds
to the raw matrix does not recreate cellranger's own cell-calling, it just
applies cell-level QC thresholds to a matrix that is mostly not cells yet.

**Detection.** Check which directory was loaded. `raw_feature_bc_matrix/`
requires its own cell-calling step (e.g. `DropletUtils::emptyDrops()`) before
anything in this pack applies.

**Right.** Use `filtered_feature_bc_matrix/` by default. If there is a
specific reason to start from the raw matrix (e.g. to recover cells cellranger's
knee-point missed), say so explicitly and perform cell-calling first.

---

## P-012 — `subset()` on a v5 object with unjoined layers

**Wrong**, on a multi-sample object straight out of `split()`:

```r
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
```

**Why it looks correct.** `subset()` runs without error and returns an
object with fewer cells, which is exactly what filtering is supposed to do.

**What actually happened.** Same root cause as `R-004`/`reference/pitfalls.md#p-004`,
applied to `subset()` specifically rather than to a `GetAssayData()` call:
on an object whose `RNA` assay has been `split()` into per-sample layers and
not yet rejoined, `subset()`'s behaviour depends on the Seurat version and
layer state, and in some cases silently filters within whichever layer is
addressable rather than consistently across every sample. The failure is
quiet — a smaller-than-expected but still plausible-looking object — not a
loud error.

**Detection.**

```r
stopifnot(sum(grepl("^counts", SeuratObject::Layers(seu[["RNA"]]))) <= 1)
```

before calling `subset()`.

**Right.** `JoinLayers()` before filtering, same as before any other
single-matrix operation:

```r
seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8)
```

`filter_cells()` in `implementation/qc_filter.R` checks for this and stops
with a clear message rather than filtering silently within one layer.

---

## P-013 — `percent.mito` computed on normalised rather than count data

**Wrong**, on an object that has already been normalised for some other
reason:

```r
seu <- Seurat::NormalizeData(seu)
seu <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mito")
```

**Why it looks correct.** `PercentageFeatureSet()` runs without error and
returns a numeric column named `percent.mito`, which looks exactly like QC
output from any other order of operations.

**What actually happened.** `PercentageFeatureSet()` sums whichever layer is
current (by default `data`) and expresses one gene set's sum as a percentage
of the total. On raw counts that percentage means "share of UMIs in this
cell." On log-normalised values it is a percentage of log-transformed
numbers — not a share of anything a reader would recognise as "percent
mitochondrial," even though the column name and the code both look identical
to the correct version.

**Detection.**

```r
stopifnot("counts do not look like raw integers - QC metrics must run before normalisation" =
            all(seu$nCount_RNA == round(seu$nCount_RNA)))
```

**Right.** Compute QC metrics before any normalisation step — this is the
order the course chapter and this pack's `add_qc_metrics()` both use. If an
object arrives already normalised for an unrelated reason, recompute metrics
explicitly from the `counts` layer rather than trusting whatever the current
default layer is.

---

## P-014 — Ribosomal percentage used as a filtering criterion

**Wrong.**

```r
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mito < 8 & percent.ribo > 5)
```

presented as if `percent.ribo` were an established QC criterion on the same
footing as `percent.mito`.

**Why it looks correct.** `percent.ribo` is computed right alongside
`percent.mito` and `percent.globin` in this pack (and in the course chapter),
in the same `VlnPlot`, and reads like one of three equally-justified QC
metrics.

**What actually happened.** Source: the course lecture for this chapter,
stated directly by the instructor: *"I've never really found any good reason
to use this percentage to filter cells... quite a lot of people do, but I'm
not entirely sure why... I'd be happy to hear an explanation."* Unlike
`percent.mito` — which has a stated biological mechanism (damaged membranes
leak cytoplasmic but not mitochondrial RNA, see `reference/citations.md#c-01`)
— there is no equivalent widely-accepted mechanism establishing a
`percent.ribo` threshold as evidence of cell quality. The course's own final
`subset()` call computes `percent.ribo` but filters only on `nFeature_RNA`
and `percent.mito` — consistent with treating it as informative, not as a
filtering criterion.

**Detection.** Any `subset()`/filter condition that includes `percent.ribo`
as a threshold.

**Right.** Compute and report `percent.ribo` — it can carry real biological
signal (`reference/citations.md#c-03`) and is worth inspecting — but do not
filter cells on it unless a specific, stated reason applies to this dataset.

---

## P-015 — QC treated as a one-shot step instead of an iterative loop

**Wrong.** Choosing thresholds once from the initial per-cell distributions,
filtering, and moving on to normalisation without ever revisiting the
decision.

**Why it looks correct.** It matches the linear structure of a tutorial:
compute metrics, look at plots, pick thresholds, filter, done.

**What actually happened.** Source: the course lecture for this chapter. The
recommended practice is to filter mildly at first — removing only the
clearest outliers — then proceed through dimensionality reduction and
clustering, and colour the resulting UMAP by `percent.mito`. If a cluster
emerges that is driven mainly by high mitochondrial content (cells that don't
share a biological identity, only a high `percent.mito`, clustering
together), that is a sign the initial filtering was too permissive, and a
reason to return to the threshold, tighten it, and re-cluster. A threshold
chosen once from the unclustered distribution and never revisited can miss
this.

**Detection.** Not statically checkable — this is a process check, not a code
pattern. Ask whether clustering was inspected for a mitochondrial-driven
cluster before treating QC as finished.

**Right.** Treat the first threshold as provisional. After an initial
clustering pass, check for and address any cluster explained mainly by
`percent.mito` before finalising the filtering step.

---

*(Further pitfalls to add as they're found.)*
