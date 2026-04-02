ctla4_gene <- "Ctla4"
melanogenesis_genes <- c("Mitf","Tyr","Tyrp1","Dct","Pmel","Mlana","Sox10","Mc1r","Kit","Gpr143")
genes_of_interest <- unique(c(ctla4_gene, melanogenesis_genes))

focus_RC <- res_RC_df[res_RC_df$gene_symbol %in% genes_of_interest, ]
focus_OC <- res_OC_df[res_OC_df$gene_symbol %in% genes_of_interest, ]
focus_RO <- res_RO_df[res_RO_df$gene_symbol %in% genes_of_interest, ]

write.csv(focus_RC, "FOCUS_CTLA4_Melanogenesis_Responder_vs_Control.csv", row.names = FALSE)
write.csv(focus_OC, "FOCUS_CTLA4_Melanogenesis_Other_vs_Control.csv", row.names = FALSE)
write.csv(focus_RO, "FOCUS_CTLA4_Melanogenesis_Responder_vs_Other.csv", row.names = FALSE)

mat <- assay(vsd)
gene_symbols_all <- sub(".*_", "", rownames(mat))
idx <- which(gene_symbols_all %in% genes_of_interest)
mat_focus <- mat[idx, , drop = FALSE]
rownames(mat_focus) <- make.unique(sub(".*_", "", rownames(mat_focus)))
mat_focus_scaled <- t(scale(t(mat_focus)))

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

png("Heatmap_CTLA4_Melanogenesis.png", width = 1200, height = 1000, res = 150)
pheatmap(mat_focus_scaled, annotation_col = annotation_col, main = "CTLA4 + Melanogenesis genes")
dev.off()

mel_idx <- which(gene_symbols_all %in% melanogenesis_genes)
mel_mat <- mat[mel_idx, , drop = FALSE]
mel_cluster_score <- colMeans(mel_mat, na.rm = TRUE)

cluster_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  melanogenesis_score = mel_cluster_score
)

write.csv(cluster_df, "Melanogenesis_cluster_score.csv", row.names = FALSE)
aggregate(melanogenesis_score ~ group, cluster_df, mean)