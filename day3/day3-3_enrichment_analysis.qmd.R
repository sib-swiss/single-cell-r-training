# Enrichment analysis

# Material

# Download the presentation manually from the following link:
# https://your-website.com/assets/pdf/202503_enrichment_analysis.pdf

# Resources:
# - MSigDB: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
# - clusterProfiler vignette: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
# - Revigo: http://revigo.irb.hr/
# - SPIA: https://bioconductor.org/packages/release/bioc/html/SPIA.html
# - GSEA original paper: https://www.pnas.org/content/102/43/15545
# - STRING for protein-protein interactions: https://string-db.org/
# - GO figure! tool: https://gitlab.com/evogenlab/GO-Figure
# - GO figure! paper: https://www.frontiersin.org/articles/10.3389/fbinf.2021.638255/full

# Load necessary data

tum_vs_norm <- readRDS("day3/tum_vs_norm_day3-2.rds")
limma_de <- readRDS("day3/limma_de_day3-2.rds")
proB <- readRDS("day3/proB_day3-2.rds")

# Load required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(msigdbdf)

# List allowed label types
AnnotationDbi::keytypes(org.Hs.eg.db)

# Select downregulated genes in tumor samples
tum_down <- subset(limma_de, limma_de$logFC < -1 & limma_de$adj.P.Val < 0.05)
tum_down_genes <- rownames(tum_down)

# Perform Gene Ontology (GO) enrichment analysis
tum_vs_norm_go <- clusterProfiler::enrichGO(
  gene = tum_down_genes,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 50
)

# View results
head(tum_vs_norm_go@result[, 1:7])

# Simplify GO results to remove redundancy
enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
head(enr_go@result[, 1:7])

# Generate an enrichment map plot
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go), showCategory = 30)

# Perform enrichment analysis with Hallmark collection
gmt <- msigdbr::msigdbr(species = "human", category = "H")

tum_vs_norm_enrich <- clusterProfiler::enricher(
  gene = tum_down_genes,
  universe = rownames(proB),
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  TERM2GENE = gmt[, c("gs_name", "gene_symbol")]
)

# View significant results
head(tum_vs_norm_enrich@result[tum_vs_norm_enrich@result$p.adjust < 0.05, 1:7])

# Clear environment
rm(list = ls())
gc()
.rs.restartR()