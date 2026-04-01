# =========================
# normalized count 저장
# =========================

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene <- rownames(norm_counts_df)

# gene 컬럼을 맨 앞으로 이동
norm_counts_df <- norm_counts_df[, c("gene", setdiff(colnames(norm_counts_df), "gene"))]

write.csv(norm_counts_df, "Normalized_counts.csv", row.names = FALSE)