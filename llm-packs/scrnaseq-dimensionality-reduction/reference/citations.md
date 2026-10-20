# Citations

| id | claim | source | verified |
|---|---|---|---|
| C-DR-01 | PCA loadings (the genes RunPCA reports per PC) are the rotation matrix coefficients, not a differential expression result or correlation score | course lecture recording, `day2/day2-2_dimensionality_reduction.qmd` video, confirmed by a live audience poll | yes (course material — see note below) |
| C-DR-02 | t-SNE preserves local neighbourhood structure but not global/inter-cluster distance; UMAP preserves more global structure but still only approximately | course lecture recording | yes (course material — see note below) |
| C-DR-03 | t-SNE perplexity / UMAP n.neighbors must be re-tuned for a subsetted (smaller) dataset | course lecture recording | yes (course material — see note below) |
| C-DR-04 | t-SNE (t-distributed stochastic neighbour embedding) method | van der Maaten & Hinton 2008, *JMLR* 9:2579-2605 | **pending** |
| C-DR-05 | UMAP (Uniform Manifold Approximation and Projection) method | McInnes, Healy & Melville 2018, arXiv:1802.03426 | **pending** |
| C-DR-06 | In scRNA-seq, the first two principal components typically capture only a small fraction of total variance, unlike bulk RNA-seq | course lecture recording (instructor's own stated experience) | yes (course material — see note below; weaker provenance than a published measurement) |

## A note on C-DR-01, C-DR-02, C-DR-03 and C-DR-06

These are sourced to the course's own lecture recording, not to
independently published literature. Present them as course
guidance/instructor experience, not as an externally peer-reviewed finding -
same convention established in the other packs' citations.md files.

## Sources deliberately not cited

- distill.pub's "How to Use t-SNE Effectively" and the UMAP author's YouTube
  talk, both linked from the course chapter. Excellent explainers, but
  informal/unversioned - point to them, don't cite them as evidence for a
  specific claim.
