# Citations

Maps each substantive claim in this pack to a source, so a model using the
pack can cite without inventing. Rule: **if a claim is not in this table,
state that it is uncited.** Do not generate a plausible-looking reference.

`verified` means a maintainer has opened the source and confirmed it supports
the specific claim made here. Entries marked `pending` were carried in from
drafting and **must be checked before this pack is released**.

| id | claim | source | verified |
|---|---|---|---|
| C-INT-01 | Integration methods (MNN, Seurat CCA/anchors) work by finding corresponding cells across batches, computing a pairwise correction vector, and applying it — not by re-estimating expression from first principles | course lecture recording, `day2/day2-3_integration.qmd` video | yes (course material, not external literature — see note below) |
| C-INT-02 | Integration-corrected output should not be used as input to differential expression; batch/technology should instead be modelled as a covariate on the original data | course lecture recording (instructor's stated recommendation, explicitly framed partly as professional opinion rather than a formal result) | yes (course material) |
| C-INT-03 | Seurat v5 layer-based CCA integration method (`IntegrateLayers`) | Hao et al. 2024, *Nat Biotechnol*, "Dictionary learning for integrative, multimodal and scalable single-cell analysis" | **pending** (same claim as `scrnaseq-qc-filtering` pack's C-05 — cross-reference rather than re-verify independently once one maintainer has checked it) |
| C-INT-04 | The original Seurat anchor-based integration method (CCA + mutual-anchor correction) | Stuart et al. 2019, *Cell*, "Comprehensive Integration of Single-Cell Data" | **pending** |
| C-INT-05 | Batch-correction method comparison/benchmarking for scRNA-seq | comparison paper linked from `advanced_exercise/02_advanced_analysis.qmd`, *Sci Rep* — https://www.nature.com/articles/s41598-021-82126-3 | **pending** |
| C-INT-06 | Mutual nearest neighbours (MNN) integration method | `batchelor`/MNN original method paper (Haghverdi et al.) — **find exact citation before release** | **pending — maintainer to locate the primary MNN paper** |

## A note on C-INT-01 and C-INT-02

These two claims are sourced to the course's own lecture recording, not to
independently published literature. They are the headline reason this pack
exists (see `reference/pitfalls.md#p-int-001`), and they are marked `yes` here
because the maintainer has directly reviewed the source (the transcript), but
they should be presented to a user as **course guidance**, not as an
externally peer-reviewed finding. `SKILL.md` reflects this distinction
explicitly.

## Sources deliberately not cited

- The Seurat integration vignette at satijalab.org. It is the correct place
  to send a human, and it changes without notice, so a pack that pins
  behaviour should not cite it as evidence (same policy as the QC pack).
- Blog posts and forum answers, however good.
