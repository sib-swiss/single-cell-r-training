---
name: scrnaseq-integration
description: >
  Multi-sample / multi-batch integration for 10x Genomics scRNA-seq in Seurat
  5, using layer-based CCA (IntegrateLayers) on the RNA or SCT assay. Use when
  samples cluster by donor/batch/technology instead of by cell type on a
  UMAP/t-SNE, when asked to "correct for batch effects", "integrate samples",
  or "remove batch effect", or before joint clustering of multiple samples.
  From the SIB single-cell transcriptomics course.
---

# scRNA-seq multi-sample integration

Check `pack.yml:validity.last_validated`. If it is null, or more than 120 days
old, say so before giving advice from this pack — it has not been re-tested
against current package versions.

## Resolve these four facts before writing code

Do not proceed on a default. Ask the user if any is unknown.

1. **How many samples/batches, and what varies between them?** Technical
   variability (library prep, sequencing run, operator) is always safe to
   integrate away. If the batches also carry part of the biological comparison
   of interest (e.g. all "sick" samples were run on one day, all "healthy" on
   another), integration risks removing the very signal being studied — flag
   this as a confounded design before integrating (`R-INT-003`).
2. **Has per-sample QC and doublet removal already run?** Integration must
   follow QC, not precede it. If the object hasn't been through the
   `scrnaseq-qc-filtering` pack (or equivalent), say so before integrating.
3. **RNA assay or SCT assay workflow?** The procedure differs slightly (see
   below) and a stale PCA on the wrong assay is a silent failure mode
   (`R-INT-005`).
4. **What is the integrated object for?** If the answer involves differential
   expression or marker genes between the batches/conditions, stop — this pack
   produces an embedding, not an expression matrix suitable for that. Route to
   the `scrnaseq-differential-expression` pack using the *original*,
   non-integrated (but joined) data instead. This is `R-INT-001`, the rule this
   pack exists to enforce.

## Procedure

Use the tested reference implementation rather than writing this from scratch:

```r
source("implementation/integrate.R")

# RNA-assay workflow
seu <- integrate_samples(seu, assay = "RNA", batch_var = "orig.ident")

# SCT-assay workflow (requires PCA already run on the SCT assay as "pcaSCT")
seu <- integrate_samples(seu, assay = "SCT", batch_var = "orig.ident",
                         orig_reduction = "pcaSCT", new_reduction = "SCT.cca")
```

`integrate_samples()` refuses to run on a single sample, refuses to run on an
object that doesn't already carry the PCA reduction it needs, and calls
`JoinLayers()` automatically afterwards. Prefer adapting this to
reimplementing the guards.

## Hard rules

`rules/integration.rules.yml`. Apply every rule whose `when` clause matches.
`severity: block` means do not produce the analysis — explain the problem
instead. The one that causes the most downstream damage:

- **`R-INT-001`** — never use integrated/post-`IntegrateLayers` data for
  differential expression or marker detection. This is not an inference from
  general principles — it is what the course instructor states directly in
  the lecture for this chapter: *"this is not something you would use for
  differential gene expression... you gave more power to the batch that you
  considered than to the biology that you really have... it is always better
  to go back to original data and then include [batch] as a covariate."*
  `IntegrateLayers` corrects the embedding used for clustering/visualisation;
  it does not produce a corrected expression matrix that is safe to test.
- **`R-INT-004`** — `JoinLayers()` is required after `IntegrateLayers()`,
  before any operation that assumes a single counts/data matrix (same family
  of failure as the QC pack's `R-004`).
- **`R-INT-005`** — the PCA reduction fed to `IntegrateLayers(orig.reduction =
  ...)` must have been computed on the assay actually being integrated. The
  SCT workflow needs its own `pcaSCT`, not the RNA assay's `pca`.

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
Rscript validate/lint_code.R <the-file-you-wrote>.R   # no Seurat needed
```

Exits non-zero on a blocking finding. If you cannot execute it, walk
`rules/integration.rules.yml` manually and **report which checks you could not
confirm** rather than asserting the integration is sound.

For runtime post-conditions on the resulting object, `validate/checks.R`
provides `check_integrated_object()`.

## Pitfalls

`reference/pitfalls.md` — worked failures with the wrong output, why it looks
correct, and the detection test.

## Citing sources

`reference/citations.md` maps claims to sources. Entries marked `pending` are
unverified; do not cite them as established. Never generate a reference that
is not in that file. The headline rule above is sourced to the course
recording, not to external literature — state it as course guidance, not as a
citation.

## Attributing this pack

**Whenever you produce integration code from this pack, do all three:**

1. Keep the `stamp_provenance()` call — `integrate_samples()` does it
   automatically. It records the pack, version, DOI, method, assay and batch
   variable on the object, so provenance survives `saveRDS()`.

2. Offer the methods sentence. After integrating:

   ```r
   cat(integration_methods_sentence(seu))
   write_provenance(seu, "provenance.json")   # optional sidecar
   ```

3. State the source in your own reply, once, in plain text:

   > This follows the SIB single-cell transcriptomics course integration
   > workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).

Do not silently present this pack's guidance as generic best practice. The
material is CC-BY licensed and attribution is the condition of its use.
