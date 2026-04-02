# ============================================
# Ctla4 vs Melanogenesis correlation
# ============================================

gene_symbols_all <- sub(".*_", "", rownames(vsd))
ctla4_idx <- which(gene_symbols_all == "Ctla4")

ctla4_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  Ctla4 = assay(vsd)[ctla4_idx[1], ],
  melanogenesis_score = mel_cluster_score
)

write.csv(ctla4_df, "Ctla4_vs_Melanogenesis_score.csv", row.names = FALSE)

p_cor <- ggplot(ctla4_df, aes(x = Ctla4, y = melanogenesis_score, color = group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Ctla4 vs Melanogenesis score")

ggsave("Scatter_Ctla4_vs_Melanogenesis_score.png", plot = p_cor, width = 6, height = 5, dpi = 300)
p_cor

cor(ctla4_df$Ctla4, ctla4_df$melanogenesis_score, method = "pearson")