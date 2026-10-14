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

Use the tested reference implementation rather than writing this from scratch:

```r
source("implementation/de_markers.R")

# Per-cluster marker discovery - fine as per-cell Wilcoxon, exploratory use
markers <- cluster_markers(seu, min_pct = 0.25, only_pos = TRUE)

# Cross-condition comparison - pseudobulk, requires >=2 replicates/condition
de <- pseudobulk_de(seu, sample_col = "orig.ident", condition_col = "type",
                    replicate_col = "patient.id",   # NULL if not genuinely paired
                    contrast = c("ETV6-RUNX1", "PBMMC"))
de_diagnostic_plots(de)   # p-value histogram + MA plot - look before trusting
```

`pseudobulk_de()` refuses to run with fewer than 2 replicates per condition
and warns loudly if `replicate_col` implies pairing that doesn't match the
number of distinct values per condition.

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

1. Keep the `stamp_provenance()` call — both `cluster_markers()` and
   `pseudobulk_de()` do it automatically.
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
