# Citations

| id | claim | source | verified |
|---|---|---|---|
| C-ANN-01 | `SingleR` correlates per-cell (or per-cluster) expression against a reference's cell-type profiles to assign labels, with a fine-tuning step to resolve closely-correlated candidate labels | Aran et al. 2019, *Nat Immunol* 20:163-172 (the `SingleR` method paper) | **pending** |
| C-ANN-02 | `SingleR` has no "rejection" option; every cell receives a label even if its true type is absent from the reference | course lecture recording, `day3/day3-1_cell_annotation.qmd` video, corroborated by the method-comparison paper linked from the course chapter | yes (course material — see note below) |
| C-ANN-03 | mRNA level and surface protein level for the same gene/marker are not always tightly correlated at the single-cell level | course lecture recording (CITE-seq CD4 example) | yes (course material — see note below) |
| C-ANN-04 | Cell-type abundance proportions estimated from scRNA-seq can be biased by capture/dissociation effects and should be validated across biological replicates and, ideally, an orthogonal method | course lecture recording | yes (course material — see note below) |
| C-ANN-05 | Review/benchmark of automated cell type annotation methods for scRNA-seq | linked from the course chapter — https://www.sciencedirect.com/science/article/pii/S2001037021000192 | **pending** |
| C-ANN-06 | CellMarker 2.0 is a manually curated database of cell-type marker genes for human and mouse, usable to ground a marker-based cluster annotation in an external, versioned source rather than memorized knowledge | Hu C, Li T, Xu Y, Zhang X, Li F, Bai J, Chen J, Jiang W, Yang K, Ou Q, Li X, Wang P, Zhang Y. "CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data." *Nucleic Acids Research* 51(D1):D870-D876, 2023. DOI: 10.1093/nar/gkac947 | yes (bibliographic record independently confirmed via Crossref) |
| C-ANN-07 | PanglaoDB is a web server and curated marker-gene database for mouse and human single-cell RNA-seq, usable as a second trusted source for the same check | Franzén O, Gan LM, Björkegren JLM. "PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data." *Database* (Oxford), 2019. DOI: 10.1093/database/baz046 | yes (bibliographic record independently confirmed via Crossref) |

## A note on C-ANN-02, C-ANN-03 and C-ANN-04

Sourced to the course's own lecture recording, not to independently
published literature (beyond the corroborating benchmark paper for
C-ANN-02, which should be checked before treating that half as externally
verified). Present as course guidance/instructor experience — same
convention as the other packs.

## A note on C-ANN-06 and C-ANN-07

"Verified" here means the bibliographic record itself (authors, title,
journal, year, DOI) was independently cross-checked against Crossref, not
that this pack's `implementation/annotate.R` has been exercised against a
live download from either database — that step still depends on each
database's current website structure, which can and does change (CellMarker
migrated its own site to a different interface during this pack's own
authoring). Verify the live download link at the time of use rather than
assuming the Procedure's example URL still resolves exactly as shown.
