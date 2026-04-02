# ============================================
# CTLA4 + melanogenesis heatmap
# ============================================

library(pheatmap)

mat <- assay(vsd)

gene_symbols_all <- sub(".*_", "", rownames(mat))
idx <- which(gene_symbols_all %in% genes_of_interest)

mat_focus <- mat[idx, , drop = FALSE]
rownames(mat_focus) <- make.unique(sub(".*_", "", rownames(mat_focus)))

mat_focus_scaled <- t(scale(t(mat_focus)))

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

png("Heatmap_CTLA4_Melanogenesis.png", width = 1200, height = 1000, res = 150)
pheatmap(
  mat_focus_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  main = "CTLA4 + Melanogenesis genes"
)
dev.off()