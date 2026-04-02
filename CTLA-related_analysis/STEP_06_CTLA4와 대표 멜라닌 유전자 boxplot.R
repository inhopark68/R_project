# ============================================
# 개별 유전자 boxplot
# ============================================

library(ggplot2)

plot_gene_box <- function(gene_symbol, vsd, meta) {
  gene_symbols <- sub(".*_", "", rownames(vsd))
  idx <- which(gene_symbols == gene_symbol)

  if (length(idx) == 0) {
    message(paste("Gene not found:", gene_symbol))
    return(NULL)
  }

  df <- data.frame(
    expression = assay(vsd)[idx[1], ],
    group = meta$group
  )

  p <- ggplot(df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    theme_bw() +
    ggtitle(paste("Expression of", gene_symbol)) +
    xlab("") +
    ylab("VST expression")

  ggsave(
    filename = paste0("Boxplot_", gene_symbol, ".png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )

  return(p)
}

plot_gene_box("Ctla4", vsd, meta)
plot_gene_box("Mitf", vsd, meta)
plot_gene_box("Tyr", vsd, meta)
plot_gene_box("Tyrp1", vsd, meta)
plot_gene_box("Dct", vsd, meta)
plot_gene_box("Pmel", vsd, meta)