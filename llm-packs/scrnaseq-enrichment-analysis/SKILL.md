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
```

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

Keep the `stamp_provenance()` call (`run_ora()` does it automatically),
offer `enrichment_methods_sentence(res)`, and state the source:

> This follows the SIB single-cell transcriptomics course enrichment
> analysis workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).
