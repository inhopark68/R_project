# =========================
# Sample distance heatmap 저장
# =========================

library(pheatmap)

vsd <- vst(dds, blind = FALSE)

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