# Citations

| id | claim | source | verified |
|---|---|---|---|
| C-ENR-01 | Over-representation analysis's background/universe must reflect the genes actually tested, not the whole genome, or the Fisher's exact test comparison is biased | the course chapter's own `GeneRatio`/`BgRatio` explanation, combined with the in-course `enrichGO()`/`enricher()` inconsistency (see `reference/pitfalls.md#p-enr-001`) | yes (course material, derived from an internal inconsistency rather than a single statement) |
| C-ENR-02 | Original Gene Set Enrichment Analysis (GSEA) method | Subramanian et al. 2005, *PNAS* 102(43):15545-15550 — https://www.pnas.org/content/102/43/15545 | yes (linked directly from the course chapter) |
| C-ENR-03 | GSEA requires a full, manually pre-ranked gene list including non-significant genes; it is not computed automatically from fold changes | course lecture recording, `day3/day3-3_enrichment_analysis.qmd` video (shared with the DE pack) | yes (course material — see note below) |
| C-ENR-04 | clusterProfiler package for ORA/GSEA of GO, KEGG and custom gene sets | Wu et al. 2021, *The Innovation* 2(3):100141 (clusterProfiler vignette linked from the course chapter) | **pending** |
| C-ENR-05 | GO terms are hierarchical and share genes extensively, producing redundant "independent" hits in a raw ORA result | course chapter's use of `clusterProfiler::simplify()`, and the linked GO-Figure! paper — https://www.frontiersin.org/articles/10.3389/fbinf.2021.638255/full | **pending** |

## A note on C-ENR-01 and C-ENR-03

C-ENR-01 is derived from an internal inconsistency in the course material
rather than a single explicit statement — treat it as course-derived
guidance, not an externally peer-reviewed finding. C-ENR-03 is sourced to
the lecture recording directly, same convention as the other packs.

## Sources deliberately not cited

- Revigo and SPIA, linked from the course chapter as tools to try, not as
  sources for a specific claim.
