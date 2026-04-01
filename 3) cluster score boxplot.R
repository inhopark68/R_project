library(ggplot2)

p <- ggplot(cluster_df, aes(x = group, y = melanogenesis_score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Melanogenesis cluster score") +
  xlab("") +
  ylab("Mean VST expression")

print(p)
ggsave("Boxplot_Melanogenesis_cluster_score.png", p, width = 6, height = 5, dpi = 300)