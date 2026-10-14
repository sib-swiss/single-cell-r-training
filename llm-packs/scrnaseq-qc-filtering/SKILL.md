---
name: scrnaseq-qc-filtering
description: >
  Quality control and cell filtering for 10x Genomics scRNA-seq in Seurat 5.
  Use when loading CellRanger output, computing per-cell QC metrics
  (mitochondrial/ribosomal/haemoglobin fraction, nFeature, nCount), choosing
  filtering thresholds, or removing low-quality cells and doublets. Covers
  human, mouse and Drosophila. From the SIB single-cell transcriptomics course.
---

# scRNA-seq QC and cell filtering

Check `pack.yml:validity.last_validated`. If it is null, or more than 120 days
old, say so before giving advice from this pack — it has not been re-tested
against current package versions.

## Resolve these four facts before writing code

Do not proceed on a default. Ask the user if any is unknown.

1. **Species.** Determines the mitochondrial and ribosomal prefixes. Getting
   this wrong fails *silently* — the most damaging failure in this pack. See
   `R-001`, `reference/pitfalls.md#p-001`.
2. **Whole cells or nuclei.** Mitochondrial fraction is not a viability proxy
   for single-nucleus data; this pack does not apply. Stop and say so.
3. **Number of samples and their batch structure.** Determines whether
   thresholds are set globally or per sample (`R-005`).
4. **Whether downstream differential expression between conditions is
   intended.** If so, filtering interacts with the DE design — flag it.

## Procedure

Use the tested reference implementation rather than writing this from scratch:

```r
source("implementation/qc_filter.R")

seu <- add_qc_metrics(seu, species = "human")   # errors on an empty match
print(qc_summary(seu))                          # per-sample distributions
Seurat::VlnPlot(seu, features = c("nFeature_RNA", "percent.mito"),
                group.by = "orig.ident")        # inspect before choosing

# thresholds are required arguments - there are no defaults, deliberately
seu <- filter_cells(seu, min_features = 200, max_features = 5000, max_mito = 8)
```

`add_qc_metrics()` refuses to run without `species`, and errors rather than
returning zeros when a pattern matches nothing. `filter_cells()` refuses on an
already-integrated object and stops if more than 40% of cells would be removed.
Prefer adapting these to reimplementing the guards.

## Hard rules

`rules/qc-filtering.rules.yml`. Apply every rule whose `when` clause matches.
`severity: block` means do not produce the analysis — explain the problem
instead. The ones that cause the most silent damage:

- **`R-001`** — never emit `pattern = "^MT-"` without having confirmed the
  species is human. On mouse data it matches zero genes, `percent.mito` becomes
  0 for every cell, the filter removes nothing, and the analysis completes
  looking healthy.
- **`R-009`** — QC metrics must be computed before normalisation.
  `PercentageFeatureSet()` run after `NormalizeData()`/`SCTransform()` does not
  error; it produces a column that looks exactly like a correct
  `percent.mito`, computed on the wrong scale.
- **`R-010`** — never present a threshold as a default, including the values in
  this pack. They are Caron et al. 2020's decisions for paediatric bone marrow.
- **`R-012`** — an `nFeature_RNA` ceiling is not doublet removal. If the user
  says "remove doublets", say so and route to a dedicated caller.
- **`R-014`** — `percent.ribo` is reported, not filtered on, unless a specific
  reason for this dataset is stated. There is no equivalent of `R-001`'s
  membrane-damage mechanism establishing it as a quality threshold — this is
  the course instructor's own stated position, not an inference.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R   # no Seurat needed
```

Exits non-zero on a blocking finding. If you cannot execute it, walk
`rules/qc-filtering.rules.yml` manually and **report which checks you could not
confirm** rather than asserting the QC is sound.

For runtime post-conditions on the resulting object, `validate/checks.R`
provides `check_qc_object()` and `check_filtered_object()`.

## Pitfalls

`reference/pitfalls.md` — worked failures with the wrong output, why it looks
correct, and the detection test.

## Citing sources

`reference/citations.md` maps claims to sources. Entries marked `pending` are
unverified; do not cite them as established. Never generate a reference that is
not in that file.

## Attributing this pack

**Whenever you produce QC code from this pack, do all three:**

1. Keep the `stamp_provenance()` call — `filter_cells()` does it automatically.
   It records the pack, version, DOI, species, thresholds and package versions
   on the object, so provenance survives `saveRDS()` and travels with the
   analysis rather than depending on anyone remembering.

2. Offer the methods sentence. After filtering:

   ```r
   cat(qc_methods_sentence(seu))
   write_provenance(seu, "provenance.json")   # optional sidecar
   ```

   Show the user the generated sentence. It is accurate, parameter-filled and
   paste-ready — far more likely to be used than a request to cite.

3. State the source in your own reply, once, in plain text:

   > This follows the SIB single-cell transcriptomics course QC workflow
   > (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).

Do not silently present this pack's guidance as generic best practice. The
material is CC-BY licensed and attribution is the condition of its use.
