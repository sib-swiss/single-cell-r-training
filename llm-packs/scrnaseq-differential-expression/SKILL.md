---
name: scrnaseq-differential-expression
description: >
  Differential gene expression for scRNA-seq in Seurat 5: per-cluster marker
  discovery (FindAllMarkers/FindMarkers, Wilcoxon) and cross-condition testing
  via pseudobulk aggregation with edgeR/limma. Use when asked to find marker
  genes, compare gene expression between conditions/groups/patients, run
  FindMarkers or FindAllMarkers, or interpret a log fold-change or p-value
  from a single-cell experiment. From the SIB single-cell transcriptomics
  course.
---

# scRNA-seq differential gene expression

Check `pack.yml:validity.last_validated`. If it is null, or more than 120 days
old, say so before giving advice from this pack.

## Resolve these facts before writing code

Do not proceed on a default. Ask the user if any is unknown.

1. **Marker discovery (per-cluster) or cross-condition comparison?** These
   need different tools and different caution. Marker discovery
   (`FindAllMarkers`) is exploratory and per-cell testing is fine for it.
   Cross-condition comparison (disease vs. healthy, treated vs. control)
   needs fact 2 before anything else.
2. **How many independent biological replicates per condition?** The single
   most important fact for cross-condition testing. Cells from the same
   patient/mouse/sample are not independent replicates of each other — see
   `R-DE-001`. If the answer is fewer than 2 per condition, say the
   comparison cannot be done with proper replication (`not_applicable_to` in
   `pack.yml`).
3. **Is this a genuinely paired/matched design, or does it only look like
   one?** Do not assume pairing because a template or prior script used it.
   See `R-DE-002` — this course's own exercise contains a worked example of
   exactly this trap.
4. **More than 2 groups?** If so, a sequence of pairwise `FindMarkers` calls
   under-corrects for multiple testing — see `R-DE-004`. Route to an
   ANOVA-like omnibus test (`limma`'s F-test) instead, or explicitly
   re-adjust the combined p-values afterwards.
5. **Publication or exploratory?** Publication-bound results need the
   diagnostic plots in `R-DE-005` before anything is reported as trustworthy.

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

Write the actual steps directly into the user's own script — the same
functions this course taught (`FindAllMarkers`, `AggregateExpression`,
`DGEList`, `voom`, `lmFit`, ...) — rather than sourcing
`implementation/de_markers.R` and calling `cluster_markers()`/
`pseudobulk_de()` as opaque wrapper functions. Those functions are the
CI-tested version of this exact logic, matching the guards already listed
below — open the file only if something fails in a way those rules don't
explain. Don't call them in the file you hand back.

Marker discovery is simple enough to stay a one-liner:

```r
markers <- Seurat::FindAllMarkers(seu, min.pct = 0.25, only.pos = TRUE, verbose = FALSE)
markers <- subset(markers, markers$p_val_adj < 0.05)
```

Cross-condition comparison is the course's own multi-step pseudobulk recipe,
written out in full rather than hidden behind a function call — long, but
every line is something this course already taught:

```r
sample_col <- "orig.ident"; condition_col <- "type"; replicate_col <- "patient.id"
contrast <- c("ETV6-RUNX1", "PBMMC")

# R-DE-001: a cross-condition comparison needs >=2 independent replicates
# per condition - this is not a per-cell test, on purpose.
md <- seu@meta.data
n_per_condition <- table(unique(md[, c(sample_col, condition_col)])[[condition_col]])
stopifnot("fewer than 2 replicates per condition" = all(n_per_condition[contrast] >= 2))

# R-DE-002: replicate_col implies a paired design - confirm with the user
# that samples are GENUINELY matched before trusting a paired result. This
# course's own exercise fabricates pairing for teaching purposes only.

bulk <- Seurat::AggregateExpression(seu, group.by = sample_col,
                                    return.seurat = TRUE, assay = "RNA", verbose = FALSE)
meta <- unique(md[, c(sample_col, condition_col, replicate_col)])
rownames(meta) <- meta[[sample_col]]
bulk@meta.data <- meta[colnames(bulk), , drop = FALSE]

counts <- as.matrix(Seurat::GetAssayData(bulk, layer = "counts"))
y <- edgeR::DGEList(counts, samples = bulk@meta.data)
y <- y[edgeR::filterByExpr(y, group = bulk@meta.data[[condition_col]]), ]

design <- stats::model.matrix(~0 + type + patient.id, data = bulk@meta.data)
colnames(design) <- make.names(sub("^type", "", colnames(design)))
contrast_mat <- limma::makeContrasts(contrasts = "ETV6.RUNX1 - PBMMC", levels = design)

v <- limma::voom(edgeR::calcNormFactors(y), design)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design), contrast_mat))
de <- limma::topTable(fit, number = Inf, sort.by = "P")
de$gene <- rownames(de)
```

(If `condition_col`/`replicate_col` have names that aren't valid R
identifiers on their own, see `implementation/de_markers.R:pseudobulk_de()`
for the more defensive column-renaming this skips for readability.)

For 3+ groups, don't repeat the pairwise block above — use an omnibus F-test
instead (`R-DE-004`): same shape, but `design <- model.matrix(~0 + group_col)`
and `limma::topTable(fit, sort.by = "F")` with no contrast. See
`implementation/de_markers.R:anova_de()` for the full pattern.

```r
source("implementation/de_markers.R")   # for diagnostics + attribution helpers only
de_diagnostic_plots(de)   # p-value histogram + MA plot - look before trusting, R-DE-005
de <- stamp_provenance(de, step = "pseudobulk_de",
                       contrast = paste(contrast, collapse = " vs "),
                       n_replicates = as.list(n_per_condition[contrast]),
                       paired = TRUE)   # FALSE if replicate_col was not used
```

## Hard rules

`rules/differential-expression.rules.yml`. `severity: block` means do not
produce the analysis — explain the problem instead. The headline:

- **`R-DE-001`** — never test a cross-condition comparison with a per-cell
  test (`FindMarkers`/`FindAllMarkers`) without flagging pseudoreplication.
  The course's own exercise is the worked example: `FindMarkers` between
  tumor and healthy pro-B cells returns 3820 "significant" genes; the
  patient-aware pseudobulk `limma` analysis on the *same* comparison returns
  far fewer. The exercise poses this as a neutral "are you surprised?"
  question. This pack exists to close that gap: route cross-condition
  comparisons to `pseudobulk_de()`, which aggregates to one value per
  biological replicate before testing.
- **`R-DE-002`** — never assume a paired/matched design without confirming
  real matching. The course's own exercise contains an explicit callout
  admitting the PBMMC/ETV6-RUNX1 "patient" pairing is fabricated "for
  demonstration purposes." Carrying that fiction into a real analysis invents
  statistical power that doesn't exist.
- **`R-DE-004`** — a sequence of pairwise `FindMarkers` calls across more than
  two groups does not control the family-wise error rate the way a single
  omnibus test does; each call's `p_val_adj` is computed only within that
  comparison. Sourced from the lecture, not written in the exercise itself.

## Write the code to a file

Write the code you produce to an actual file before running or checking it —
see [`../ENVIRONMENT.md`](../ENVIRONMENT.md#write-the-code-to-a-file-always)
for why this applies to every pack and is not optional.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R   # no Seurat needed
```

For runtime post-conditions, `validate/checks.R` provides
`check_pseudobulk_design()` and `check_de_diagnostics()`.

## Pitfalls

`reference/pitfalls.md` — worked failures, including the course's own
pseudoreplication numbers as the concrete illustration.

## Citing sources

`reference/citations.md` maps claims to sources. Never generate a reference
that is not in that file.

## Attributing this pack

**Whenever you produce DE code from this pack, do all three:**

1. Call `stamp_provenance()` explicitly, as shown above, passing the
   `step`/`contrast`/`n_replicates`/`paired` you actually used (or
   `step = "cluster_markers"` for marker discovery).
2. Offer the methods sentence after running DE:

   ```r
   cat(de_methods_sentence(de))
   write_provenance(de, "provenance.json")
   ```

3. State the source in your own reply, once, in plain text:

   > This follows the SIB single-cell transcriptomics course differential
   > expression workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).

Do not silently present this pack's guidance as generic best practice. The
material is CC-BY licensed and attribution is the condition of its use.
