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

Candidates for packs 2 and 3, in order of pitfall density: integration
(`day2-3`) and differential expression (`day3-2`).

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

## Running the checks

Static linting needs no Seurat and runs anywhere:

```bash
Rscript llm-packs/scrnaseq-qc-filtering/validate/test-lint.R
Rscript llm-packs/scrnaseq-qc-filtering/validate/lint_code.R scripts/run_analysis.R
```

Full validation needs the course container:

```bash
docker run --rm -v "$PWD":/work -w /work \
  sibswiss/training-singlecell-rstudio:2026.6 \
  Rscript llm-packs/scrnaseq-qc-filtering/validate/run_all.R
```

CI runs both weekly and on changes under `llm-packs/` — see
`.github/workflows/validate-llm-packs.yml`.

## Maintenance rule

Every pack has one named maintainer in `pack.yml`, not a team alias. A pack
whose validation has been red for two weeks is unpublished rather than fixed
later. Better three packs green than eight stale.
