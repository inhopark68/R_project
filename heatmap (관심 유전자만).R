library(pheatmap)

mat <- assay(vsd)

# 관심 유전자 index
gene_idx <- rownames(mat)[sub(".*_", "", rownames(mat)) %in% genes_of_interest]

mat_goi <- mat[gene_idx, , drop = FALSE]

# gene symbol로 row 이름 변경
rownames(mat_goi) <- sub(".*_", "", rownames(mat_goi))

# scale
mat_goi_scaled <- t(scale(t(mat_goi)))

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

pheatmap(
  mat_goi_scaled,
  annotation_col = annotation_col,
  main = "CTLA4 + Melanogenesis genes"
)