# ============================================
# Melanogenesis cluster score boxplot
# ============================================

p_cluster <- ggplot(cluster_df, aes(x = group, y = melanogenesis_score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Melanogenesis cluster score") +
  xlab("") +
  ylab("Mean VST expression")

ggsave("Boxplot_Melanogenesis_cluster_score.png", plot = p_cluster, width = 6, height = 5, dpi = 300)
p_cluster