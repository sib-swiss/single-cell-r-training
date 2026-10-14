# Pitfalls

---

## P-DR-001 — PCA loadings misread as differential expression or correlation

**Wrong.** "`RunPCA`'s output shows which genes are differentially expressed
along PC1" or "...which genes correlate with PC1."

**Why it looks correct.** `RunPCA` prints a list of genes with positive and
negative values for each PC, and genes with strong, consistent differences
between groups *also* tend to show up as strong loadings — the two concepts
correlate in practice often enough that conflating them doesn't immediately
produce an obviously wrong plot.

**What actually happened.** Source: the course lecture, directly confirmed
by a live audience poll during the session — asked whether the genes
`RunPCA` lists per PC are (a) differentially expressed genes, (b)
correlation scores, or (c) the rotation matrix, a majority of participants
initially answered (a) or (b). Both are wrong. A PC is, by construction, a
linear combination of the original genes: `PC1 = a1*Gene1 + a2*Gene2 + ...`.
The values `RunPCA` reports are the `a` coefficients (the loadings/rotation
matrix) — how much each gene contributes to that axis — not a hypothesis
test result and not a correlation coefficient. No p-value, no test, no
comparison group is involved.

**Detection.** Any explanation of PCA loadings that uses the words
"differentially expressed," "significant," or "correlated with" the PC.

**Right.** "Genes listed by `RunPCA` for a PC are the ones contributing most
strongly, positively or negatively, to that axis's definition (the rotation
matrix) — not a statistical comparison between groups of cells."

---

## P-DR-002 — Skipping PCA before t-SNE/UMAP

**Wrong.** Running `RunUMAP()`/`RunTSNE()` on a raw or merely-scaled
expression matrix (e.g. via a custom embedding path, or on data imported
from outside Seurat) without a PCA step first.

**Why it looks correct.** Both algorithms accept any numeric matrix as
input; nothing in their signature requires PCA specifically.

**What actually happened.** Source: course lecture, stated directly: PCA
removes redundancy (correlated genes collapse onto the same axis) and keeps
only the directions carrying real variance, before t-SNE/UMAP ever see the
data. Skipping it means t-SNE/UMAP optimise over thousands of raw, highly
correlated and noisy gene dimensions instead of ~25-50 informative ones —
slower, and the resulting layout is not the well-behaved "clusters now
visually separable" result the course's exercises show.

**Detection.**

```r
stopifnot("pca" %in% names(seu@reductions))
```

before `RunUMAP()`/`RunTSNE()`.

**Right.** `RunPCA()` first, then `RunUMAP(seu, dims = 1:N)` — always
specifying `dims` from the PCA, never running UMAP on the full gene set.

---

## P-DR-003 — Reading cluster distance off a t-SNE plot

**Wrong.** "Cluster A and cluster B are close together in the t-SNE plot, so
they must be transcriptionally similar" (or the reverse: "far apart, so very
different").

**Why it looks correct.** It's a 2D scatterplot; reading spatial distance as
similarity is the natural way to read any scatterplot, and it is the correct
reading for PCA.

**What actually happened.** Source: course lecture. t-SNE's optimisation
target only preserves *local* neighbourhood relationships (who is a
near-neighbour of whom) — global, between-cluster distance is not part of
what it optimises for, and is not recoverable from the plot. The course's
own demonstration makes this concrete: given a toy dataset where two of
three clusters are genuinely closer together than the third, t-SNE plots all
three as roughly equidistant regardless of perplexity, while UMAP does
partially recover the true relative distances (though the course notes this
recovery is inconsistent) — but only once the clusters are far enough apart
in sample count for the effect to show at all in their demonstration.

**Detection.** Any claim about relative similarity between (not within)
clusters based on a t-SNE plot's layout.

**Right.** Use distance *within* a cluster (local structure) as the thing
t-SNE plots show meaningfully. For an actual quantitative comparison between
clusters, use marker genes / differential expression, not plot geometry. If
using UMAP, note that global distance is only approximately preserved, not
exact.

---

## P-DR-004 — Default neighbourhood parameter reused on a subsetted dataset

**Wrong.** Subsetting to one cell type (e.g. only T cells) and running
`RunUMAP(seu_subset, dims = 1:25)` with the default `n.neighbors = 30`
unchanged.

**Why it looks correct.** The function call is identical to the one used on
the full dataset and runs without error or warning.

**What actually happened.** Source: course lecture: "if one has subset the
data... the number of neighbors in a cluster would anyway be most of the
time lower than 30, then 30 is too much." A neighbourhood size tuned for a
dataset of thousands of cells across many cell types does not automatically
make sense for a subset with far fewer cells per expected substructure.

**Detection.** `n.neighbors`/`perplexity` left at the default after a
`subset()` call earlier in the same analysis.

**Right.** Lower the neighbourhood parameter for a subsetted analysis, and
say why.

---

## P-DR-005 — A specific `dims=` range stated as correct with no elbow plot

**Wrong.** "Use `dims = 1:25` for the UMAP" with no supporting plot or
reasoning, for a dataset this specific choice wasn't derived from.

**Why it looks correct.** `dims = 1:25` is literally the value used in this
course's own exercises, so it looks like received wisdom.

**What actually happened.** It is this course's example dataset's own elbow
plot result (Caron et al. 2020 bone marrow). The course's own "Answer C"
exercise shows this directly: `dims = 1:5` gives an undifferentiated "blob,"
`dims = 1:50` over-separates, and `dims = 1:100` cannot even run because
`RunPCA` only computed 50 PCs by default — the right number depends on the
elbow plot (or another diagnostic) for the dataset actually being analysed,
not a transplanted constant.

**Detection.** A `dims=` value used with no elbow plot (or equivalent)
generated for *this* object.

**Right.** Generate `ElbowPlot(seu, ndims = 40)` (or more) for the object at
hand, and choose `dims` based on where it flattens — and if more precision
is needed than the default 50 PCs allow, increase `npcs` at the `RunPCA`
step, not `dims` past what was computed.
