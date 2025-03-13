library(dplyr)

# get synonyms from flybase
# https://flybase-ftp.s3.us-east-1.amazonaws.com/releases/FB2025_01/precomputed_files/genes/fbgn_annotation_ID_fb_2025_01.tsv.gz
syn <- read.delim("fbgn_annotation_ID_fb_2025_01.tsv", skip = 4) |>
  rename(gene_symbol = X..gene_symbol) |>
  mutate(annotation_ID = paste("Dmel", annotation_ID, sep = "_")) |>
  select(gene_symbol, annotation_ID)

datadirs <- list.files(path = "data", full.names = TRUE)
feat_files <- file.path(datadirs, "features.tsv.gz")

#feat_file <- "data/Female_Cocaine_1/features.tsv.gz"

for (feat_file in feat_files) {
  
  # join synonym file with feature file
  feat <- gzfile(feat_file) |>
    read.delim(header = FALSE) |>
    rename(gene_id = V1,
           gene_id2 = V2,
           type = V3) |>
    left_join(syn, by = join_by(gene_id == annotation_ID)) |>
    select(gene_id, gene_symbol, type) |>
    mutate(gene_symbol = ifelse(is.na(gene_symbol),
                                gene_id,
                                gene_symbol))
  
  # # remove old file
  # file.remove(feat_file)
  
  # write new feature file
  gz_out <- gzfile(feat_file, "w")
  write.table(
    feat,
    gz_out,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
  close(gz_out)
  
}

sparse_matrix <- Seurat::Read10X(data.dir = datadirs[1])

seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "CocaineStudy")
