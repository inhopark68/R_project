# MA plot, sample distance heatmap, normalized counts 저장
# gene symbol 반영 heatmap도 포함

# ============================================
# STEP 5_추가시각화_05.R
# MA plot, sample distance heatmap, normalized counts 저장
# ============================================

library(DESeq2)
library(pheatmap)

# STEP 4에서 사용한 annotation DB와 동일하게 설정
library(AnnotationDbi)

species_db <- "mouse"

if (species_db == "mouse") {
  library(org.Mm.eg.db)
  anno_db <- org.Mm.eg.db
} else {
  library(org.Hs.eg.db)
  anno_db <- org.Hs.eg.db
}

# 결과 다시 준비
res_RC <- results(dds, contrast = c("group", "Responder", "Control"))
res_OC <- results(dds, contrast = c("group", "Other", "Control"))
res_RO <- results(dds, contrast = c("group", "Responder", "Other"))

# VST 변환
vsd <- vst(dds, blind = FALSE)

# ============================================
# 1) MA plot
# ============================================

png("MAplot_Responder_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(res_RC, main = "MA plot: Responder vs Control", ylim = c(-5, 5))
dev.off()

png("MAplot_Other_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(res_OC, main = "MA plot: Other vs Control", ylim = c(-5, 5))
dev.off()

png("MAplot_Responder_vs_Other.png", width = 1200, height = 1000, res = 150)
plotMA(res_RO, main = "MA plot: Responder vs Other", ylim = c(-5, 5))
dev.off()

# ============================================
# 2) Sample distance heatmap
# ============================================

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

sample_labels <- paste(colnames(vsd), colData(vsd)$group, sep = "_")
rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels

png("Sample_distance_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample-to-sample distances"
)
dev.off()

# ============================================
# 3) Normalized counts 저장
# ============================================

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene_id <- rownames(norm_counts_df)
norm_counts_df$gene_id_clean <- sub("\\..*$", "", norm_counts_df$gene_id)

norm_counts_df$gene_symbol <- mapIds(
  anno_db,
  keys = norm_counts_df$gene_id_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

norm_counts_df$gene_name <- mapIds(
  anno_db,
  keys = norm_counts_df$gene_id_clean,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

norm_counts_df <- norm_counts_df[, c(
  "gene_id",
  "gene_id_clean",
  "gene_symbol",
  "gene_name",
  setdiff(colnames(norm_counts_df), c("gene_id", "gene_id_clean", "gene_symbol", "gene_name"))
)]

write.csv(norm_counts_df, "Normalized_counts.csv", row.names = FALSE)

# ============================================
# 4) Top variable genes heatmap
# gene symbol 반영
# ============================================

mat <- assay(vsd)
gene_var <- apply(mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, nrow(mat))]
mat_top <- mat[top_genes, , drop = FALSE]

top_genes_clean <- sub("\\..*$", "", rownames(mat_top))
top_symbols <- mapIds(
  anno_db,
  keys = top_genes_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

row_labels <- ifelse(is.na(top_symbols) | top_symbols == "", rownames(mat_top), top_symbols)
rownames(mat_top) <- make.unique(row_labels)

mat_top_scaled <- t(scale(t(mat_top)))

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

png("Heatmap_top50_variable_genes.png", width = 1200, height = 1000, res = 150)
pheatmap(
  mat_top_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  main = "Top 50 variable genes"
)
dev.off()

# ============================================
# 5) 완료 메시지
# ============================================

cat("\n===== STEP 5 완료 =====\n")
cat("저장된 파일:\n")
cat("- MAplot_Responder_vs_Control.png\n")
cat("- MAplot_Other_vs_Control.png\n")
cat("- MAplot_Responder_vs_Other.png\n")
cat("- Sample_distance_heatmap.png\n")
cat("- Normalized_counts.csv\n")
cat("- Heatmap_top50_variable_genes.png\n")