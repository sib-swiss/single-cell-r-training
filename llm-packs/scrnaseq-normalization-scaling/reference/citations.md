# Citations

| id | claim | source | verified |
|---|---|---|---|
| C-NORM-01 | LogNormalize (per-cell total-count scaling + log transform) accounts for sequencing-depth variation between cells | course chapter `day2-1_normalization_scaling.qmd` | yes (course material) |
| C-NORM-02 | SCTransform uses a more sophisticated normalization model and is argued to perform better than LogNormalize, at the cost of speed and transparency | Hafemeister & Satija 2019, *Genome Biology* 20:296 — https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1 | yes (cited in the course) |
| C-NORM-03 | `ScaleData()` defaults to scaling only `VariableFeatures(object)`, not every gene in the assay | Seurat documentation (`?Seurat::ScaleData`) | **pending** — verify against the installed Seurat version before marking verified; this is version-sensitive API behaviour, not a citable scientific claim |

## Sources deliberately not cited

- The Seurat vignettes at satijalab.org, for the same reason as the other
  packs: unversioned, changes without notice.
