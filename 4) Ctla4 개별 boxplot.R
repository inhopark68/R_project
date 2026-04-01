ctla4_idx <- which(sub(".*_", "", rownames(vsd)) == "Ctla4")

df_ctla4 <- data.frame(
  expression = assay(vsd)[ctla4_idx[1], ],
  group = meta$group
)

p_ctla4 <- ggplot(df_ctla4, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Ctla4 expression") +
  xlab("") +
  ylab("VST expression")

print(p_ctla4)
ggsave("Boxplot_Ctla4.png", p_ctla4, width = 6, height = 5, dpi = 300)