---
name: scrnaseq-enrichment-analysis
description: >
  Gene set enrichment analysis for scRNA-seq DE results: over-representation
  analysis (ORA, clusterProfiler::enrichGO/enricher) and GSEA. Use when
  running enrichment analysis, interpreting GO/KEGG/MSigDB/Hallmark results,
  choosing an ORA universe/background, or preparing a ranked gene list for
  GSEA. From the SIB single-cell transcriptomics course.
---

# scRNA-seq enrichment analysis

Check `pack.yml:validity.last_validated`. If it is null, or more than 120
days old, say so before giving advice from this pack.

## Resolve these facts before writing code

1. **ORA (gene list only) or GSEA (full ranked list with statistics)?**
   GSEA needs every tested gene's statistic, not just the significant ones —
   see `R-ENR-003`.
2. **What is the correct background/universe for ORA?** The genes actually
   detectable/tested in this experiment, not the whole genome — see
   `R-ENR-001`. This is the headline rule: the course's own two ORA calls in
   the same chapter disagree on this.
3. **Up-regulated, down-regulated, or both, tested separately?** Pooling both
   directions into one undirected gene list loses the direction of effect —
   see `R-ENR-002`.
4. **Which gene ID type?** `clusterProfiler` needs `keyType` to match what the
   gene list actually uses (symbol, Ensembl, Entrez).

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
tum_down_genes <- rownames(subset(limma_de, logFC < -1 & adj.P.Val < 0.05))

go_res <- clusterProfiler::enrichGO(
  gene          = tum_down_genes,
  universe      = rownames(proB),       # R-ENR-001: not the clusterProfiler default
  OrgDb         = "org.Hs.eg.db",
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 50)
go_res <- clusterProfiler::simplify(go_res)   # R-ENR-006: reduce GO term redundancy

source("implementation/enrichment.R")   # for the attribution helper only
go_res <- stamp_provenance(go_res, n_genes = length(tum_down_genes),
                           n_universe = length(rownames(proB)), ontology = "BP")
```

This pack also has a tested `run_ora()` wrapper in `implementation/enrichment.R`
that refuses to run without an explicit `universe`, matching `R-ENR-001`
above — open it only if something fails in a way that rule doesn't explain.
Write the steps out as shown, not a call to it.

## Hard rules

`rules/enrichment-analysis.rules.yml`. The headline rule, directly evidenced
by an inconsistency in the course's own material:

- **`R-ENR-001`** — ORA's `universe=` must be the genes actually
  tested/detectable in the experiment. The course's own `enrichGO()` call
  for this chapter omits `universe=` (so clusterProfiler defaults to every
  annotated gene in `org.Hs.eg.db`), while its own later `enricher()` call in
  the *same* chapter correctly sets `universe = rownames(proB)`. The chapter
  explains exactly why this matters (`GeneRatio`/`BgRatio`) without applying
  it consistently in its own first example.
- **`R-ENR-003`** — GSEA requires the full ranked gene list, every gene with
  a statistic, ranked *manually* beforehand — `GSEA()`-family functions do
  not rank it for you.

## Write the code to a file

Write the code you produce to an actual file before running or checking it —
see [`../ENVIRONMENT.md`](../ENVIRONMENT.md#write-the-code-to-a-file-always)
for why this applies to every pack and is not optional.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R
```

`validate/checks.R` provides `check_enrichment_result()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Call `stamp_provenance()` explicitly as shown above, then offer the methods
sentence and append it to the project's running methods file:

```r
cat(enrichment_methods_sentence(res))
append_methods(enrichment_methods_sentence(res))
```

State the source:

> This follows the SIB single-cell transcriptomics course enrichment
> analysis workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
