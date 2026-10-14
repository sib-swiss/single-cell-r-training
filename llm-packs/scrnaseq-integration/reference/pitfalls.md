# Pitfalls

Each entry: the wrong output, **why it looks correct**, what actually happened,
how to detect it, the right version. Written as worked failures rather than
advice, following the convention established in the
`scrnaseq-qc-filtering` pack.

---

## P-INT-001 — Integrated data used for differential expression

**Wrong.**

```r
seu <- Seurat::IntegrateLayers(seu, method = CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca")
seu[["RNA"]] <- Seurat::JoinLayers(seu[["RNA"]])
Seurat::DefaultAssay(seu) <- "integrated"   # or similar
markers <- Seurat::FindMarkers(seu, ident.1 = "disease", ident.2 = "healthy")
```

**Why it looks correct.** The object integrated cleanly, `FindMarkers` runs
without error, and a lot of pre-Seurat-v5 tutorials online do exactly this —
older workflows build a corrected `integrated` assay and use it for
everything downstream, including DE. The code is idiomatic for *that* older
pattern.

**What actually happened.** `IntegrateLayers` computes a correction vector in
a reduced space so that corresponding cells from different batches sit closer
together for clustering/visualisation. It does not produce expression values
that preserve the magnitude of real biological differences. Testing on
integrated data gives the batch correction more statistical weight than the
biology — the course instructor's own words, from the lecture recording for
this chapter: *"these are quite harsh, the numbers that you get at the end
are not really representative anymore of what happened at the beginning."*
The result can understate a real effect, or manufacture one, depending on how
the batch structure relates to the comparison being tested.

**Detection.** Any DE/marker call whose input object or assay was produced by
`IntegrateLayers` (directly, or via `DefaultAssay<-` to an assay derived from
it) is a finding, regardless of whether it also "looks reasonable."

**Right.**

```r
seu <- Seurat::IntegrateLayers(seu, method = CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca")
seu[["RNA"]] <- Seurat::JoinLayers(seu[["RNA"]])
Seurat::DefaultAssay(seu) <- "RNA"          # original, joined, non-integrated

# simple case: Wilcoxon on original data, batch noted as a caveat
markers <- Seurat::FindMarkers(seu, ident.1 = "disease", ident.2 = "healthy")

# complex case: batch/technology modelled explicitly as a covariate
# (route to the scrnaseq-differential-expression pack's pseudobulk workflow)
```

---

## P-INT-002 — Treating a confounded design as a pure technical batch effect

**Wrong.** Integrating all "disease" samples (collected in one run) against
all "healthy" samples (collected in another run), then reporting that
integration "worked" because the UMAP now overlaps nicely, without checking
whether the disease signal survived.

**Why it looks correct.** The before/after UMAP comparison is the course's own
headline visual for "why integrate" — batches that used to separate now
overlap, which is the expected, desired outcome when the only difference
between samples truly is technical.

**What actually happened.** When the batch structure is confounded with the
biological comparison of interest, the integration method cannot tell "same
cell type, different batch, please merge" from "genuinely different cell
state caused by the condition, please don't merge." A UMAP that now overlaps
nicely can mean the correction removed noise — or that it removed signal.

**Detection.** Before integrating, ask whether batches/samples differ only
technically, or also carry part of the biological comparison. After
integrating, check a positive control: a gene or cell-type proportion known
(from the literature, or from the pre-integration data) to differ between
conditions should still show that difference post-integration. If it has
flattened out, say so rather than reporting a clean UMAP as success.

**Right.** State the design (confounded or not) explicitly before
integrating, and verify a positive control survives integration before
treating the result as trustworthy for downstream interpretation (not for DE —
see P-INT-001 — but for cluster composition / annotation purposes).

---

## P-INT-003 — Stale or wrong-assay PCA fed into `IntegrateLayers`

**Wrong**, when integrating the SCT assay:

```r
seu[["SCT_v5"]] <- as(seu[["SCT"]], Class = "Assay5")
seu[["SCT_v5"]] <- split(seu[["SCT_v5"]], f = seu$orig.ident)
seu <- Seurat::IntegrateLayers(seu, assay = "SCT_v5", method = CCAIntegration,
                               orig.reduction = "pca",   # <- RNA-assay PCA
                               new.reduction = "SCT.cca")
```

**Why it looks correct.** `seu[["pca"]]` exists, the function call has no
missing arguments, and it runs to completion.

**What actually happened.** `pca` was computed on the RNA assay's variable
features, not the SCT assay's residuals. The anchors `IntegrateLayers`
computes are only as good as the space they're computed in; reusing the wrong
PCA silently integrates on a basis that does not reflect the data actually
being corrected.

**Detection.** Before integrating assay `X`, confirm the reduction named in
`orig.reduction` was computed with `RunPCA` on assay `X`, after its current
variable features were selected — not reused from a different assay or an
earlier run.

**Right.**

```r
Seurat::DefaultAssay(seu) <- "SCT"
seu <- Seurat::RunPCA(seu, reduction.name = "pcaSCT", verbose = FALSE)
seu[["SCT_v5"]] <- as(seu[["SCT"]], Class = "Assay5")
seu[["SCT_v5"]] <- split(seu[["SCT_v5"]], f = seu$orig.ident)
seu <- Seurat::IntegrateLayers(seu, assay = "SCT_v5", method = CCAIntegration,
                               orig.reduction = "pcaSCT",
                               new.reduction = "SCT.cca")
```

---

*(Further pitfalls to add as they're found: integrating before re-running
variable feature selection per layer; mixing MNN/Harmony embeddings with
CCA-integrated ones in the same downstream analysis; label-transfer used as a
substitute for proper integration when the reference and query differ
substantially in composition.)*
