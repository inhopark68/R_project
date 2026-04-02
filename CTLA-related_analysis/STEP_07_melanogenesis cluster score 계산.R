# ============================================
# Melanogenesis cluster score
# ============================================

mat <- assay(vsd)
gene_symbols_all <- sub(".*_", "", rownames(mat))

mel_idx <- which(gene_symbols_all %in% melanogenesis_genes)
mel_mat <- mat[mel_idx, , drop = FALSE]
rownames(mel_mat) <- sub(".*_", "", rownames(mel_mat))

mel_cluster_score <- colMeans(mel_mat, na.rm = TRUE)

cluster_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  melanogenesis_score = mel_cluster_score
)

write.csv(cluster_df, "Melanogenesis_cluster_score.csv", row.names = FALSE)

aggregate(melanogenesis_score ~ group, cluster_df, mean)