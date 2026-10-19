---
name: scrnaseq-cell-annotation
description: >
  Cell type annotation for scRNA-seq in Seurat 5: manual marker-gene
  inspection and automated reference-based annotation with SingleR/celldex.
  Use when annotating clusters, identifying cell types, running SingleR,
  choosing a reference, or interpreting marker gene / module score results.
  From the SIB single-cell transcriptomics course.
---

# scRNA-seq cell type annotation

Check `pack.yml:validity.last_validated`. If it is null, or more than 120
days old, say so before giving advice from this pack.

## Resolve these facts before writing code

1. **Manual, automated, or both?** Automated needs a reference that actually
   contains the expected cell types. Manual needs marker genes for each
   candidate cell type fetched from a trusted database (CellMarker,
   PanglaoDB, or a connected tool) — never assigned from the model's own
   memorized knowledge as the default path (`R-ANN-007`). Course lecture's
   own answer: try both and compare, when feasible.
2. **Does a suitable reference exist for this organism/tissue?** If not
   (non-model organism, unusual tissue), automated annotation will still
   assign *a* label to every cell — it has no "none of the above" option
   with `SingleR`. See `R-ANN-002`.
3. **Is the question about cell-type ABUNDANCE differing between
   conditions?** That needs multiple biological replicates, not a single
   sample's proportions — see `R-ANN-005`.
4. **Which assay/layer is being correlated against the reference?** `slot=`
   is deprecated in Seurat 5 — see `R-ANN-001`.

## Procedure

> **Missing a package?** See [`../ENVIRONMENT.md`](../ENVIRONMENT.md) before
> installing anything ad hoc - don't resolve missing dependencies one at a
> time as you discover them.

```r
ref <- celldex::NovershternHematopoieticData()     # or another suitable reference
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu),  # default layer="data", not slot=
                                ref = ref, labels = ref$label.main)

# R-ANN-004: relabel rare (<10-cell) labels to "none", stated explicitly,
# rather than leaving them cluttering every plot.
tab <- table(seu_SingleR$labels)
rare <- names(tab)[tab < 10]
final_labels <- seu_SingleR$labels
final_labels[final_labels %in% rare] <- "none"
seu$SingleR_annot <- final_labels

SingleR::plotScoreHeatmap(seu_SingleR)              # inspect before trusting labels - R-ANN-002
SingleR::plotDeltaDistribution(seu_SingleR)

source("implementation/annotate.R")   # for the attribution helper only
seu <- stamp_provenance(seu, n_labels = length(unique(seu_SingleR$labels)),
                        n_rare_relabelled = sum(seu_SingleR$labels %in% rare))
```

This pack also has a tested `annotate_with_reference()` wrapper in
`implementation/annotate.R` doing the same thing, matching the guards
already listed above — open it only if something fails in a way those don't
explain. Write the steps out as shown, not a call to it.

**Manual, marker-based annotation** — ground every cell-type call in a
fetched marker table, never in memorized/recalled gene-to-cell-type
knowledge (`R-ANN-007`). Only fall back to recalled knowledge if no trusted
source is reachable, and if you do, say so explicitly in your reply, not
only in a code comment:

```r
markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, verbose = FALSE)
top_markers <- do.call(rbind, by(markers, markers$cluster,
                                 function(d) head(d[order(d$p_val_adj), ], 10)))

# Fetch a trusted, versioned marker reference. Prefer a connected tool if
# your assistant has one (e.g. a ToolUniverse/MCP CellMarker lookup) - it
# stays current with no download link to track. Otherwise download the
# database's own current table directly (verify the live link on the page
# first - these restructure occasionally):
#   CellMarker: https://bio-bigdata.hrbmu.edu.cn/CellMarker/ (Excel download)
#   PanglaoDB:  https://panglaodb.se/markers.html            (tab-separated)
marker_db <- readxl::read_excel(here::here("Cell_marker_Human.xlsx"))
colnames(marker_db)           # inspect before assuming a schema - it has changed before
celltype_col <- "cell_name"   # set from colnames() above, not guessed
gene_col     <- "Symbol"      # set from colnames() above, not guessed

gene_sets <- split(marker_db[[gene_col]], marker_db[[celltype_col]])
for (cl in unique(top_markers$cluster)) {
  cl_genes <- top_markers$gene[top_markers$cluster == cl]
  overlap <- vapply(gene_sets, function(g) length(intersect(cl_genes, g)), integer(1))
  cat("Cluster", cl, "- best database match:", names(sort(overlap, decreasing = TRUE))[1],
      "(", max(overlap), "of", length(cl_genes), "top markers overlap)\n")
}

Seurat::FeaturePlot(seu, features = top_markers$gene[top_markers$cluster == levels(Seurat::Idents(seu))[1]])

source("implementation/annotate.R")   # for match_clusters_to_markers() + the attribution helper only
seu <- stamp_provenance(seu, marker_source = "CellMarker")   # or "PanglaoDB" / "connected_tool" /
                                                              # "model_memory_unverified" if that was the fallback
```

The overlap above is a starting point to confirm against `FeaturePlot`, not a
label to accept unexamined — and `implementation/annotate.R` has a tested
`match_clusters_to_markers()` version of the same computation, taking the
marker table's column names as explicit arguments precisely because their
names aren't something to hardcode blindly.

## Hard rules

`rules/cell-annotation.rules.yml`. The one with a concrete known instance in
this repository:

- **`R-ANN-001`** — `slot=` is deprecated in Seurat 5. **Known instance:**
  `scripts/run_analysis.R:209,296` still calls `SingleR::SingleR(test =
  Seurat::GetAssayData(seu_int, slot = "data"), ...)`, while the rendered
  `day3-1_cell_annotation.qmd` chapter correctly uses the layer-implicit
  default. The script and the chapter disagree.
- **`R-ANN-002`** — every cell gets a `SingleR` label, even if its true type
  isn't in the reference. Say so before trusting a full annotation.
- **`R-ANN-005`** — cell-type abundance differences between conditions need
  multiple biological replicates before being reported as real, not a
  single-sample proportion comparison.
- **`R-ANN-007`** — a marker-based cell-type call must be grounded in
  CellMarker, PanglaoDB, or a connected tool backed by one of them — never
  presented as evidence-based when it was actually just recalled from
  training knowledge. Memory is a last resort only, used explicitly and
  flagged as such, never the default path.

## Write the code to a file

Write the code you produce to an actual file before running or checking it —
see [`../ENVIRONMENT.md`](../ENVIRONMENT.md#write-the-code-to-a-file-always)
for why this applies to every pack and is not optional.

## Check your own output

```bash
Rscript validate/lint_code.R <the-file-you-wrote>.R
```

`validate/checks.R` provides `check_annotation_object()`.

## Pitfalls

`reference/pitfalls.md`.

## Citing sources

`reference/citations.md`. Never generate a reference not in that file.

## Attributing this pack

Call `stamp_provenance()` explicitly as shown above, offer
`annotation_methods_sentence(seu)`, and state the source:

> This follows the SIB single-cell transcriptomics course cell annotation
> workflow (DOI: 10.5281/zenodo.5703106, CC-BY-4.0).

For marker-based annotation specifically, also state which marker source the
labels rest on — `annotation_methods_sentence()` reflects `marker_source`
when it was passed to `stamp_provenance()`. If it was
`"model_memory_unverified"`, that sentence alone is not enough: say so in
plain text in your reply too, per `R-ANN-007`.
