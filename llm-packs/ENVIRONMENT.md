# Environment and shared conventions

Every pack's `SKILL.md` points here for things that apply identically to all
of them: what to do about a missing package, the one rule that applies to
every pack's Procedure section, how to run the resulting script without
flooding the conversation with its console output, and how to keep the
result reproducible on a different machine or a different day.

## Missing a package?

Every pack's `pack.yml:depends_on` names the packages its advice depends on.
If you hit a missing-package error while following any pack, **don't install
packages one at a time as you discover them** — that produces a separate,
often very long installation log (compiler output, every transitive
dependency's own download) for every single package, repeated once per
pack as the analysis moves through the pipeline. On a cold environment with
none of this pre-installed, that can mean the same flood six or seven times
over in one session.

### Do this instead

1. Check whether this project already has an `renv.lock` (the course repo's
   own, or a copy fetched alongside the packs). If so, restore it **once**,
   before continuing with the analysis:

   ```r
   install.packages("renv")
   renv::restore(lockfile = "renv.lock")
   # or, if only the course's own lockfile is available and not this project's:
   renv::restore(lockfile = "https://raw.githubusercontent.com/sib-swiss/single-cell-r-training/refs/heads/master/renv.lock")
   ```

   This installs every package every pack depends on in one step, at the
   exact versions each `pack.yml:validity.validated_against` was tested
   against — not a version discovered ad hoc at install time.

2. If `renv::restore()` isn't available or appropriate (e.g. no lockfile
   reachable, or only a single package is actually needed), install only
   what's missing, and batch it into **one** call covering everything
   `pack.yml:depends_on.in_lockfile` across the packs in play names, rather
   than one call per package as each gap is discovered:

   ```r
   BiocManager::install(c("Seurat", "SingleR", "celldex", "clusterProfiler", ...),
                        update = FALSE, ask = FALSE)
   ```

3. Either way, don't let the raw installation log flow back into the
   conversation unexamined — redirect it (e.g. to a file, or capture and
   only surface the last few lines / any error) and report success or
   failure, not the full trying-URL-by-URL transcript.

4. A package in `pack.yml:depends_on.recommended_but_absent` is not part of
   the course's `renv.lock` by design (see that pack's `pack.yml` for why) —
   don't add it via step 1 or 2; that's a separate, deliberate decision to
   flag to the user, not a gap to silently fill.

## Write the code to a file, always

Before anything else in any pack's Procedure: write the code you produce to
an actual file in the user's project (e.g. `analysis/<step>.R`, or appended
to whatever script they're already building) — don't just execute it in a
transient/ephemeral session with nothing saved. This is not optional. It is
what makes every pack's `validate/lint_code.R` check possible at all, and
it's what keeps the analysis inspectable, version-controlled, and
re-runnable after this conversation ends.

## Run scripts without flooding context back in

Running the code from any pack's Procedure produces console output
(package messages, printed diagnostics, warnings) that gets read back in
full once the run finishes - this is a separate cost from the installation
logs covered above, and it recurs every time a script runs, not just once
per package. Redirect it and pull in only what's needed, rather than
letting an entire run flow back unexamined:

```bash
Rscript analysis/02_dimred.R > analysis/02_dimred.log 2>&1
tail -n 20 analysis/02_dimred.log        # or: grep -i error analysis/02_dimred.log
```

If the run succeeded, a short tail (or "no errors found") is enough to
report back - the full log stays on disk if it's needed for debugging
later, rather than needing to be re-run just to produce it again. If it
failed, show the actual error from the log, not just "it failed."

## Reproducible paths, seeds, and sessions

A script that only runs correctly from one specific working directory, or
only gives the same numbers on one specific run, is not something the user
can hand to a colleague or re-run next month and trust. `here` is already in
this course's `renv.lock` - use it.

1. **Build every file path with `here::here()`, never a hardcoded absolute
   path, and never `setwd()`.** `here::here("data", "counts")`,
   `here::here("results", "qc_umap.png")` - resolved from the project root
   regardless of where the script is actually run from. A path that is
   correct only because the working directory happens to match whoever
   wrote it is not reproducible, and `setwd()` inside a script changes state
   for everything that runs after it in the same session.

2. **Set a seed once, near the top of the script, as defense-in-depth.**
   Seurat's own key functions already default to a fixed internal seed
   (`RunPCA`/`RunUMAP`: `seed.use = 42`; `RunTSNE`: `seed.use = 1`;
   `FindClusters`: `random.seed = 0`) - so results from those specific calls
   are already deterministic run-to-run on the same package versions,
   without any extra code. A top-of-script `set.seed(1)` doesn't change
   that, but it does cover anything else in the script that doesn't carry
   its own seed argument (a `sample()` call, a bootstrap, a quick manual
   subsample) - set one regardless of whether the pack's own calls strictly
   require it.

3. **Record the session alongside the output**, not only on the R object
   itself. Every pack's `stamp_provenance()` already records the R and
   package versions used *on the object* - but that's lost if only a
   derived plot or table is kept. Alongside any saved output:

   ```r
   writeLines(capture.output(sessionInfo()), here::here("results", "sessionInfo.txt"))
   ```

4. **Every pack's `append_methods()` writes to the same shared file**,
   `here::here("METHODS.md")` by default - not a per-pack file. Across a
   project that moves through several packs (QC, then normalization, then
   dimensionality reduction, ...), this builds up into one running methods
   section covering every step actually run, each keyed so re-running a
   step updates its own entry instead of duplicating it. It is a companion
   to `provenance.json`, not a replacement: `provenance.json` is the
   per-step machine-readable record on one object; `METHODS.md` is the
   cumulative, human/manuscript-facing one across the whole project.

5. **Don't use interactive-only calls** (`menu()`, `readline()`, a bare
   `View()`) in a script meant to be re-run unattended later - they block
   on input that won't be there.
