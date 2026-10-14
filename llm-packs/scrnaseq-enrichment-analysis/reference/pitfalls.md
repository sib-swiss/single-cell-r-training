# Pitfalls

---

## P-ENR-001 — ORA run against the whole-genome default background

**Wrong.**

```r
tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
                                            "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)
```

with no `universe=` argument — which is exactly the course's own exercise
code for this chapter.

**Why it looks correct.** It runs without error or warning, returns a
plausible-looking enrichment table, and is literally what the written
course chapter shows as the worked example.

**What actually happened.** With no `universe=`, `enrichGO()` defaults the
background to every gene annotated in `org.Hs.eg.db` — the whole genome's
GO annotations, not the genes this experiment could actually have detected.
The same chapter explains, in its own callout box, exactly what `GeneRatio`
and `BgRatio` mean and why the comparison matters — and then its own *later*
call in the same chapter, `enricher()` against the Hallmark collection,
correctly sets `universe = rownames(proB)`. The two calls in one chapter
disagree on this. Single-cell RNA-seq only detects a subset of the genome
per cell type (dropout), so testing against the full genome as background
systematically distorts which gene sets look over- or under-represented.

**Detection.** Any `enrichGO()`/`enricher()`/`enrichKEGG()` call with no
`universe=` argument.

**Right.**

```r
tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
                                            "org.Hs.eg.db",
                                            universe = rownames(proB),  # genes actually tested
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)
```

---

## P-ENR-002 — Up- and down-regulated genes pooled into one undirected test

**Wrong.** Combining all significant genes (both directions) from a DE
result into one gene list and running a single ORA on it.

**Why it looks correct.** It's one gene list, one function call, one result
table — simpler to produce and to read than two.

**What actually happened.** `scripts/run_analysis.R` in this repository runs
`enrichGO()` separately for `tum_down_genes` and `tum_up_genes`; the written
chapter only demonstrates the down-regulated direction as its single worked
example, which risks a learner generalising "just test the significant gene
list" without separating direction. A pathway enriched among a pooled
up+down gene list could be driven entirely by one direction, the other, or
a genuine mix — the pooled result can't distinguish these, and "this
pathway is enriched for differential expression" is a much weaker,
harder-to-interpret finding than "this pathway is enriched among
up-regulated genes."

**Right.** Run ORA on up-regulated and down-regulated gene sets separately
(as `scripts/run_analysis.R` does), and report direction alongside any
enriched term.

---

## P-ENR-003 — GSEA run on a gene list instead of a ranked statistic

**Wrong.**

```r
clusterProfiler::GSEA(geneList = significant_genes_only, ...)
```

or passing an unranked vector of all genes.

**Why it looks correct.** `GSEA()`'s first argument is also called
"geneList," which reads like it wants the same kind of gene list `enrichGO()`
takes.

**What actually happened.** Source: course lecture, stated directly: "the
GSEA function will not rank your fold changes automatically... you have to
do that manually before you use the function." GSEA's enrichment score is
computed by walking down a list ranked by a per-gene statistic (fold change
or test statistic) for **every tested gene**, not just the significant ones
— it needs the position of non-significant genes in the ranking to compute
where the enrichment score peaks. A gene list of significant genes only, or
an unranked vector, does not have the information GSEA's algorithm needs.

**Detection.** A `GSEA()`-family call where the input wasn't explicitly
sorted by a statistic beforehand, or was filtered to significant genes only.

**Right.**

```r
ranked <- sort(setNames(limma_de$logFC, rownames(limma_de)), decreasing = TRUE)
clusterProfiler::GSEA(geneList = ranked, ...)
```

---

## P-ENR-004 — A long raw GO hit list reported as independent findings

**Wrong.** Reporting 80 significantly enriched GO terms as 80 separate
discoveries, with no mention that they overlap.

**Why it looks correct.** Each term has its own row, its own adjusted
p-value, and technically each one did pass the significance threshold
independently.

**What actually happened.** Gene Ontology terms are hierarchical and
share genes extensively — many "significant" terms are near-duplicates of
each other (parent/child terms, or terms describing the same underlying
process from different angles). The course's own exercise calls
`clusterProfiler::simplify()` specifically to collapse this redundancy
before looking at results, and the lecture separately discusses `REVIGO`
and `GO-Figure!` as dedicated tools for the same purpose.

**Right.** Call `clusterProfiler::simplify()` (or an equivalent grouping
step) before reporting a GO-based ORA/GSEA result, and say that raw term
counts overstate the number of independent findings.
