# Citations

| id | claim | source | verified |
|---|---|---|---|
| C-ANN-01 | `SingleR` correlates per-cell (or per-cluster) expression against a reference's cell-type profiles to assign labels, with a fine-tuning step to resolve closely-correlated candidate labels | Aran et al. 2019, *Nat Immunol* 20:163-172 (the `SingleR` method paper) | **pending** |
| C-ANN-02 | `SingleR` has no "rejection" option; every cell receives a label even if its true type is absent from the reference | course lecture recording, `day3/day3-1_cell_annotation.qmd` video, corroborated by the method-comparison paper linked from the course chapter | yes (course material — see note below) |
| C-ANN-03 | mRNA level and surface protein level for the same gene/marker are not always tightly correlated at the single-cell level | course lecture recording (CITE-seq CD4 example) | yes (course material — see note below) |
| C-ANN-04 | Cell-type abundance proportions estimated from scRNA-seq can be biased by capture/dissociation effects and should be validated across biological replicates and, ideally, an orthogonal method | course lecture recording | yes (course material — see note below) |
| C-ANN-05 | Review/benchmark of automated cell type annotation methods for scRNA-seq | linked from the course chapter — https://www.sciencedirect.com/science/article/pii/S2001037021000192 | **pending** |

## A note on C-ANN-02, C-ANN-03 and C-ANN-04

Sourced to the course's own lecture recording, not to independently
published literature (beyond the corroborating benchmark paper for
C-ANN-02, which should be checked before treating that half as externally
verified). Present as course guidance/instructor experience — same
convention as the other packs.
