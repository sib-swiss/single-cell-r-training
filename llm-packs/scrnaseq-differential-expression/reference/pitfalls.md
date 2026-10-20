# Pitfalls

Each entry: the wrong output, **why it looks correct**, what actually happened,
how to detect it, the right version. Written as worked failures rather than
advice, following the convention established in the `scrnaseq-qc-filtering`
pack.

---

## P-DE-001 — Pseudoreplication: per-cell testing across patients

**Wrong.**

```r
seu <- Seurat::SetIdent(seu, value = "type")   # "ETV6-RUNX1" vs "PBMMC"
de <- Seurat::FindMarkers(seu, ident.1 = "ETV6-RUNX1", ident.2 = "PBMMC")
de_sig <- subset(de, de$p_val_adj < 0.05)
```

**Why it looks correct.** `FindMarkers` runs without error, returns a long,
plausible-looking table with fold changes and adjusted p-values, and this is
*literally the code in the course's own exercise* (`day3-2
differential_gene_expression.qmd`), presented as one of two valid approaches
to compare.

**What actually happened.** Every cell is treated as an independent
observation, but cells from the same patient are correlated with each other
far more than cells from different patients — patient identity itself drives
a lot of expression variation. Running a test as if there were hundreds of
independent replicates per condition, when there are really only a handful of
patients, inflates significance. The course's own numbers make this
concrete: `FindMarkers` on this exact tumor-vs-healthy pro-B cell comparison
returns **3820** significant genes; a patient-aware pseudobulk `limma`
analysis of the identical comparison returns far fewer. The exercise asks
"are you surprised by this number?" without stating which number should be
trusted — this pitfall is that gap.

**Why this looks like a reasonable thing to do:** the lecture's own analogy
makes the structure clear: if you want to compare species counts between
forest and grassland ecosystems, sampling many quadrats within one forest
site and treating each quadrat as a replicate overstates your true sample
size — the forest *site* is the actual independent replicate, not each
square metre sampled within it. A patient is the "site"; cells are the
quadrats.

**Detection.** Any cross-condition comparison test that operates at the
per-cell level, where the conditions being compared span multiple
samples/patients, without aggregation to one value per
replicate-per-condition first.

**Right.**

```r
bulk <- Seurat::AggregateExpression(seu, group.by = "sample",
                                    return.seurat = TRUE, assay = "RNA")
# one row per gene, one column per (patient x condition) - not per cell
y <- edgeR::DGEList(Seurat::GetAssayData(bulk, layer = "counts"),
                    samples = bulk@meta.data)
y <- y[edgeR::filterByExpr(y, group = bulk$type), ]
design <- model.matrix(~0 + bulk$type)
v <- limma::voom(edgeR::calcNormFactors(y), design)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design),
  limma::makeContrasts(contrasts = "...", levels = design)))
```

---

## P-DE-002 — Treating fabricated pairing as real

**Wrong.**

```r
proB$patient.id <- gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id <- sapply(strsplit(proB$patient.id, "-"), `[`, 2)
design <- model.matrix(~0 + proB$type + proB$patient.id)   # paired design
```

copied onto a *new* dataset, without checking whether its samples are
genuinely patient-matched.

**Why it looks correct.** This is exactly the course's own exercise code, it
runs, and a paired design is statistically more powerful than an unpaired
one, so reusing the pattern "looks like" good practice.

**What actually happened.** The course's own exercise contains an explicit
callout admitting this: the healthy (PBMMC) and malignant (ETV6-RUNX1)
samples were **not** patient-matched in reality — the exercise pairs
PBMMC-1 with ETV6-RUNX1-1 "as if" they came from the same patient purely "for
demonstration purposes," to show how a paired design is coded. A model that
reuses this exact code pattern on a real dataset, or that silently assumes
any `patient.id`-like column implies real pairing, invents statistical power
the data doesn't have — the paired model's smaller standard errors are not
justified.

**Detection.** Before including any blocking/pairing covariate, confirm with
the user that the pairing reflects real experimental design (the same
biological subject measured under both conditions), not a convenient
same-index assumption.

**Right.** Ask. If pairing isn't real, drop the blocking covariate and use
the unpaired pseudobulk design; say explicitly that this is less powerful but
not fabricated.

---

## P-DE-003 — Repeated pairwise `FindMarkers` across more than two groups

**Wrong.**

```r
de_young_old   <- Seurat::FindMarkers(seu, ident.1 = "young",  ident.2 = "old")
de_young_mid   <- Seurat::FindMarkers(seu, ident.1 = "young",  ident.2 = "middle")
de_mid_old     <- Seurat::FindMarkers(seu, ident.1 = "middle", ident.2 = "old")
# report each comparison's own p_val_adj as "significant at FDR 5%"
```

**Why it looks correct.** Each individual call is a perfectly valid pairwise
test, and each returns its own properly adjusted p-values — adjusted *within
that one comparison*.

**What actually happened.** Source: course lecture, not written anywhere in
the exercise itself. Running three separate pairwise tests and reporting each
one's own FDR control is not the same as controlling the error rate across
all three comparisons considered together — the instructor states this
directly when discussing a 3-group (newborn/middle-age/old mice) design: an
F-test/ANOVA-style omnibus test (available via `limma`) handles this
correctly in one pass; doing it via repeated pairwise `FindMarkers` calls
does not, and "you would have to further manually adjust for the multiple
comparisons you're doing."

**Detection.** More than one `FindMarkers` call used to cover all pairs of a
3+ level grouping variable, with no subsequent correction across the combined
set of p-values.

**Right.** Use `limma`'s F-test for the omnibus question ("does this gene
differ among any of the three groups?"), then pairwise contrasts /
post-hoc tests only as a follow-up, with p-values adjusted across the full
set of comparisons performed.

---

## P-DE-004 — Fold-change sign reported without stating direction

**Wrong.** "CD4 has a negative log2FC" with no further context.

**Why it looks correct.** It's a true statement about the number in the
output table.

**What actually happened.** `FindMarkers(ident.1 = X, ident.2 = Y)` defines
the sign relative to `ident.1`. The course's own CD8-vs-CD4 worked example
shows this is easy to misread: a negative log2FC for the CD4 gene means lower
expression in CD8+ T cells (`ident.1`) than CD4+ T cells (`ident.2`) — not
"downregulated" in some absolute sense.

**Detection.** Any reported fold-change sign with no accompanying statement
of which group was `ident.1` vs. `ident.2`.

**Right.** Always pair a fold-change sign with "higher/lower in `<ident.1
label>` relative to `<ident.2 label>`."
