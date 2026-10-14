# Environment: what to do about a missing package

Every pack's `pack.yml:depends_on` names the packages its advice depends on.
If you hit a missing-package error while following any pack, **don't install
packages one at a time as you discover them** — that produces a separate,
often very long installation log (compiler output, every transitive
dependency's own download) for every single package, repeated once per
pack as the analysis moves through the pipeline. On a cold environment with
none of this pre-installed, that can mean the same flood six or seven times
over in one session.

## Do this instead

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
