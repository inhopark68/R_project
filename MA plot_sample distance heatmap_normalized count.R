library(pheatmap)

# VST 변환
vsd <- vst(dds, blind = FALSE)

# -----------------
# MA plots
# -----------------
png("MAplot_Responder_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(
  results(dds, contrast = c("group", "Responder", "Control")),
  main = "MA plot: Responder vs Control",
  ylim = c(-5, 5)
)
dev.off()

png("MAplot_Other_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(
  results(dds, contrast = c("group", "Other", "Control")),
  main = "MA plot: Other vs Control",
  ylim = c(-5, 5)
)
dev.off()

png("MAplot_Responder_vs_Other.png", width = 1200, height = 1000, res = 150)
plotMA(
  results(dds, contrast = c("group", "Responder", "Other")),
  main = "MA plot: Responder vs Other",
  ylim = c(-5, 5)
)
dev.off()

# -----------------
# Sample distance heatmap
# -----------------
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- paste(colnames(vsd), colData(vsd)$group, sep = "_")
colnames(sample_dist_matrix) <- paste(colnames(vsd), colData(vsd)$group, sep = "_")

png("Sample_distance_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample-to-sample distances"
)
dev.off()

# -----------------
# Normalized counts 저장
# -----------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene <- rownames(norm_counts_df)
norm_counts_df <- norm_counts_df[, c("gene", setdiff(colnames(norm_counts_df), "gene"))]

write.csv(norm_counts_df, "Normalized_counts.csv", row.names = FALSE)

cat("저장 완료\n")
cat("- MAplot_Responder_vs_Control.png\n")
cat("- MAplot_Other_vs_Control.png\n")
cat("- MAplot_Responder_vs_Other.png\n")
cat("- Sample_distance_heatmap.png\n")
cat("- Normalized_counts.csv\n")