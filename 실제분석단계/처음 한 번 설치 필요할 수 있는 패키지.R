install.packages("pheatmap")
BiocManager::install(c("GEOquery", "DESeq2", "AnnotationDbi", "org.Mm.eg.db"))

사람 데이터로 바꿔야 하면:

BiocManager::install("org.Hs.eg.db")

그리고 species_db <- "human"으로 바꾸면 됩니다.