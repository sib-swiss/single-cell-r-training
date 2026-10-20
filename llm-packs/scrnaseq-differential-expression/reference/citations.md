# Citations

Maps each substantive claim in this pack to a source. Rule: **if a claim is
not in this table, state that it is uncited.** Do not generate a
plausible-looking reference.

`verified` means a maintainer has opened the source and confirmed it supports
the specific claim made here. Entries marked `pending` must be checked before
this pack is released.

| id | claim | source | verified |
|---|---|---|---|
| C-DE-01 | Per-cell differential expression tests between conditions spanning multiple samples/patients suffer inflated false-positive rates due to pseudoreplication; pseudobulk aggregation with bulk RNA-seq tools (edgeR/limma) gives better-calibrated results | Soneson & Robinson 2018, *Nat Methods* 15:255-261, "Bias, robustness and scalability in single-cell differential expression analysis" — https://www.nature.com/articles/nmeth.4612 (linked directly from `day3-2_differential_gene_expression.qmd` and described in detail in the course lecture as "an evaluation of 36 methods... still valid") | **pending** (directly linked and discussed by the course; maintainer should still open and confirm the specific claim before marking verified) |
| C-DE-02 | limma-voom ("roomLima" in the lecture transcript — likely a mishearing of "voom-limma") ranks highly among DE methods evaluated in Soneson & Robinson 2018, including for scalability | Soneson & Robinson 2018 (as C-DE-01) | **pending** |
| C-DE-03 | A healthy p-value histogram from a well-calibrated DE test is approximately L-shaped (spike near 0, flat elsewhere) | linked from the course exercise: "How to interpret a p-value histogram" — http://varianceexplained.org/statistics/interpreting-pvalue-histogram/ | **pending — this is a blog post, consider finding a peer-reviewed equivalent or explicitly citing it as an informal source, not a journal citation** |
| C-DE-04 | `muscat` as a purpose-built package for pseudobulk DE across multiple cell types and samples | linked from the course exercise's Material section — https://bioconductor.org/packages/release/bioc/html/muscat.html | **pending** |

## Sources deliberately not cited

- The OSCA book chapter on multi-sample comparisons
  (https://bioconductor.org/books/3.17/OSCA.multisample/multi-sample-comparisons.html),
  linked from the course exercise. It is the correct place to send a human for
  a fuller treatment, but as a living online book it can change without
  notice — point to it, don't pin specific claims to it.
- The ecology sampling-design analogy (forest vs. grassland quadrats) used in
  the lecture to explain pseudoreplication. It is a teaching device, not a
  citable source — the actual statistical claim it illustrates is C-DE-01.
