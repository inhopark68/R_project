melanogenesis_genes <- c(
  "Mitf", "Tyr", "Tyrp1", "Dct",
  "Pmel", "Mlana", "Sox10", "Mc1r", "Kit", "Gpr143"
)

mat <- assay(vsd)
gene_symbols_all <- sub(".*_", "", rownames(mat))

mel_idx <- which(gene_symbols_all %in% melanogenesis_genes)
mel_mat <- mat[mel_idx, , drop = FALSE]

melanogenesis_score <- colMeans(mel_mat, na.rm = TRUE)

cluster_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  melanogenesis_score = melanogenesis_score
)

aggregate(melanogenesis_score ~ group, cluster_df, mean)