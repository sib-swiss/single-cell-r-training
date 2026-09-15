# Citations

Maps each substantive claim in this pack to a source, so a model using the pack
can cite without inventing. Rule: **if a claim is not in this table, state that
it is uncited.** Do not generate a plausible-looking reference.

`verified` means a maintainer has opened the source and confirmed it supports
the specific claim made here. Entries marked `pending` were carried in from
drafting and **must be checked before this pack is released** — an unverified
citation in a pack whose whole purpose is trustworthiness is worse than none.

| id | claim | source | verified |
|---|---|---|---|
| C-01 | A damaged cell membrane loses cytoplasmic RNA faster than mitochondrial RNA, so a high mitochondrial fraction indicates damaged cells | Lun, McCarthy & Marioni 2016, *F1000Research* 5:2122 — https://f1000research.com/articles/5-2122 | yes (cited in `day1-3_analysis_tools_qc.qmd`) |
| C-02 | The threshold values used in this course (<8% mitochondrial counts, >200 features/cell) and the PBMMC/B-ALL dataset itself | Caron et al. 2020, *Sci Rep* 10:8079 — https://www.nature.com/articles/s41598-020-64929-x | yes (cited in the course; thresholds in its M&M) |
| C-03 | Ribosomal protein transcript abundance can carry biological signal rather than only technical variation | Caron et al. 2020 (as C-02) | yes |
| C-04 | MALAT1 is highly abundant in 10x poly-A data and may have biological relevance, so removal is a decision not a default | Shaat et al. 2021, *Cell Death Discov* — https://www.nature.com/articles/s41420-020-00383-y | yes (cited in the course) |
| C-05 | Seurat v5 layer-based assay structure; `layer=` replaces `slot=` | Hao et al. 2024, *Nat Biotechnol*, "Dictionary learning for integrative, multimodal and scalable single-cell analysis" | **pending** |
| C-06 | Count-based filtering does not detect homotypic doublets; dedicated callers are required | scDblFinder — Germain et al., *F1000Research* | **pending** |
| C-07 | Cardiomyocytes and other high-OXPHOS cell types legitimately exceed conventional mitochondrial thresholds | — | **pending — find a primary source or drop the specific claim** |
| C-08 | Ambient RNA contamination and its correction | Young & Behjati, SoupX, *GigaScience* 2020 | **pending** |
| C-09 | Data-driven / adaptive QC thresholds outperform fixed cutoffs | — | **pending — likely Ilicic et al. 2016, *Genome Biol*; verify it supports this** |

## Sources deliberately not cited

- The Seurat vignettes at satijalab.org. They are the correct place to send a
  human, and the course links them, but they are unversioned and change without
  notice, so a pack that pins behaviour should not cite them as evidence.
- Blog posts and forum answers, however good. If a claim only has that kind of
  support, say so — that is itself useful information for the user.
