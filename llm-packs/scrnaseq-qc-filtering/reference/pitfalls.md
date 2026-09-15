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
same problem recurs.

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

*(P-007 … P-014 follow the same shape and are still to be written: ambient RNA
not corrected; QC run after integration; global thresholds across heterogeneous
samples; MALAT1 removed without justification; raw vs filtered CellRanger matrix
confusion; `subset()` on a v5 object with unjoined layers; `percent.mito`
computed on normalised rather than count data.)*
