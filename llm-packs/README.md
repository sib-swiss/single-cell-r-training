# LLM packs

Machine-readable companions to the course chapters: the decisions, guards and
failure modes of each workflow, in a form an AI coding assistant can act on.

The course material teaches a human, linearly, with narrative and exercises. A
pack targets a different reader — a model writing analysis code for a scientist
— and so keeps the decisions and discards the pedagogy. The `.qmd` chapters stay
the human artifact; packs are generated alongside, not instead.

## Status

| Pack | Source chapter | Status |
|---|---|---|
| [`scrnaseq-qc-filtering`](scrnaseq-qc-filtering) | `day1/day1-3_analysis_tools_qc.qmd` | pilot, unreleased |
| [`scrnaseq-normalization-scaling`](scrnaseq-normalization-scaling) | `day2/day2-1_normalization_scaling.qmd` | pilot, unreleased |
| [`scrnaseq-dimensionality-reduction`](scrnaseq-dimensionality-reduction) | `day2/day2-2_dimensionality_reduction.qmd` | pilot, unreleased |
| [`scrnaseq-integration`](scrnaseq-integration) | `day2/day2-3_integration.qmd` | pilot, unreleased |
| [`scrnaseq-clustering`](scrnaseq-clustering) | `day2/day2-4_clustering.qmd` | pilot, unreleased |
| [`scrnaseq-cell-annotation`](scrnaseq-cell-annotation) | `day3/day3-1_cell_annotation.qmd` | pilot, unreleased |
| [`scrnaseq-differential-expression`](scrnaseq-differential-expression) | `day3/day3-2_differential_gene_expression.qmd` | pilot, unreleased |
| [`scrnaseq-enrichment-analysis`](scrnaseq-enrichment-analysis) | `day3/day3-3_enrichment_analysis.qmd` | pilot, unreleased |

Every course analysis chapter now has a pack except trajectory/advanced
analysis (`advanced_exercise/`), deliberately left out for now - it's
taught as a bonus overview (`eval: false`) plus one hands-on vignette,
lower pitfall density than the rest.

Course lecture recordings are mined for guidance not written in the slides or
exercises (see `pack.yml:provenance.derived_from` for each pack). Raw
auto-generated captions for all 10 course videos live under
[`_dev/transcripts/`](_dev/transcripts) as authoring material - not part of
any pack, not validated by CI.

## What a pack is

```
scrnaseq-qc-filtering/
├── pack.yml           provenance, scope limits, and the validity contract
├── SKILL.md           entry point: what to establish before writing code
├── rules/             assertions, readable by a model AND executable as tests
├── reference/         pitfalls as worked failures, and a claim→source table
├── implementation/    the CI-tested reference implementation
├── validate/          linter, runtime post-conditions, fixtures, entrypoint
└── eval/              graded tasks, run with and without the pack
```

Three design commitments are worth stating, because they are what distinguish
this from a folder of Markdown:

**A pack expires.** `pack.yml:validity` carries `last_validated`,
`validated_against` (R, Seurat and container versions) and `stale_after_days`.
Only CI writes `last_validated`, and only on a green run. A stale pack carrying
an institutional stamp is worse than no pack, because it launders a wrong answer
with our authority.

**Rules are executable.** Everything in `rules/` has a `lint` or
`post_condition` where it can, so the same file instructs a model and tests the
output. Prose warnings cannot do both.

**Claims are benefit, not assumption.** `eval/tasks.yml` measures the pack with
and without, across model versions. Until those numbers exist, "this makes
AI-assisted analysis more trustworthy" is a hypothesis.

## Writing Procedure code

Every Seurat (or edgeR/limma) call in a `SKILL.md` Procedure block or an
`implementation/*.R` reference is run, and its console output read back, by
whatever agent is following the pack to confirm each step worked. Pass
`verbose = FALSE` on every call that has it (check with
`"verbose" %in% names(formals(fn))` if unsure) — the alternative is real,
observed cost: an un-silenced `FindClusters()` swept across 8 resolutions
printed ~70 lines of Louvain/modularity output for a single step. A handful
of functions used in these packs (`RunTSNE`, `SCTransform`) have no `verbose`
argument at all; nothing to do there.

## Running the checks

The same commands apply to every pack - substitute the pack directory.
Static linting needs no Seurat and runs anywhere:

```bash
Rscript llm-packs/<pack>/validate/test-lint.R
Rscript llm-packs/<pack>/validate/lint_code.R scripts/run_analysis.R
```

Full validation needs the course container:

```bash
docker run --rm -v "$PWD":/work -w /work \
  sibswiss/training-singlecell-rstudio:2026.6 \
  Rscript llm-packs/<pack>/validate/run_all.R
```

CI discovers every pack under `llm-packs/` automatically (anything with a
`pack.yml`) and runs this for each one, weekly and on changes under
`llm-packs/` — see `.github/workflows/validate-llm-packs.yml`.

## Maintenance rule

Every pack has one named maintainer in `pack.yml`, not a team alias. A pack
whose validation has been red for two weeks is unpublished rather than fixed
later. Better three packs green than eight stale.
