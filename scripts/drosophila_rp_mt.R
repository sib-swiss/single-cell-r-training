library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Drosophila melanogaster", "EnsDb"))
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb, return.type = "data.frame") 
annotations <- annotations %>% dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
View(annotations)
mt <- annotations %>% dplyr::filter(seq_name == "mitochondrion_genome") %>% dplyr::pull(gene_name)
rp <- annotations[grepl("ribo", annotations$description, ignore.case = TRUE), ]
# mitochondrial genes:
# "^mt:"
# ribosomal genes:
# "^Rp[SL]"