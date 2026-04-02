# =========================================================
# Figure 1~5 plotting code
# =========================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(reshape2)

dir.create("Figures", showWarnings = FALSE)

# ---------------------------------------------------------
# 공통 설정
# ---------------------------------------------------------
meta$group <- factor(meta$group, levels = c("Control", "Other", "Responder"))

melanogenesis_genes <- c(
  "Mitf", "Tyr", "Tyrp1", "Dct",
  "Pmel", "Mlana", "Sox10", "Mc1r", "Kit", "Gpr143"
)

immune_genes <- c(
  "Cd8a", "Cd8b1", "Gzmb", "Prf1", "Ifng",
  "Pdcd1", "Lag3", "Tigit", "Havcr2",
  "Foxp3", "Ctla4", "Tbx21", "Stat1",
  "Cxcl9", "Cxcl10", "Cxcl11"
)

gene_symbols_all <- sub(".*_", "", rownames(vsd))
mat_vsd <- assay(vsd)

# =========================================================
# Figure 1
# Global transcriptomic landscape
# =========================================================

# 1A. group count barplot
df_group <- as.data.frame(table(meta$group))
colnames(df_group) <- c("group", "count")

p1A <- ggplot(df_group, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  ggtitle("Figure 1A. Sample groups") +
  xlab("") +
  ylab("Number of samples")

ggsave("Figures/Figure1A_group_counts.png", p1A, width = 5, height = 4, dpi = 300)

# 1B. PCA
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p1B <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("Figure 1B. PCA")

ggsave("Figures/Figure1B_PCA.png", p1B, width = 6, height = 5, dpi = 300)

# 1C. sample distance heatmap
sample_dists <- dist(t(mat_vsd))
sample_dist_matrix <- as.matrix(sample_dists)
sample_labels <- paste(colnames(vsd), colData(vsd)$group, sep = "_")
rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels

png("Figures/Figure1C_sample_distance_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Figure 1C. Sample-to-sample distances"
)
dev.off()

# 1D. DEG count summary
deg_counts <- data.frame(
  comparison = c("Responder vs Control", "Other vs Control", "Responder vs Other"),
  n_sig = c(
    sum(!is.na(res_RC$padj) & res_RC$padj < 0.05 & abs(res_RC$log2FoldChange) >= 1),
    sum(!is.na(res_OC$padj) & res_OC$padj < 0.05 & abs(res_OC$log2FoldChange) >= 1),
    sum(!is.na(res_RO$padj) & res_RO$padj < 0.05 & abs(res_RO$log2FoldChange) >= 1)
  )
)

p1D <- ggplot(deg_counts, aes(x = comparison, y = n_sig, fill = comparison)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  ggtitle("Figure 1D. Number of significant DEGs") +
  xlab("") +
  ylab("DEG count") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("Figures/Figure1D_DEG_counts.png", p1D, width = 7, height = 5, dpi = 300)

# =========================================================
# Figure 2
# CTLA4 and immune activation
# =========================================================

# 2A. Ctla4 boxplot
ctla4_idx <- which(gene_symbols_all == "Ctla4")
df_ctla4 <- data.frame(
  expression = mat_vsd[ctla4_idx[1], ],
  group = meta$group
)

p2A <- ggplot(df_ctla4, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Figure 2A. Ctla4 expression") +
  xlab("") +
  ylab("VST expression")

ggsave("Figures/Figure2A_Ctla4_boxplot.png", p2A, width = 5, height = 4.5, dpi = 300)

# 2B. immune heatmap
immune_idx <- which(gene_symbols_all %in% immune_genes)
mat_immune <- mat_vsd[immune_idx, , drop = FALSE]
rownames(mat_immune) <- make.unique(sub(".*_", "", rownames(mat_immune)))
mat_immune_scaled <- t(scale(t(mat_immune)))
mat_immune_scaled <- mat_immune_scaled[complete.cases(mat_immune_scaled), , drop = FALSE]

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

png("Figures/Figure2B_immune_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  mat_immune_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Figure 2B. Immune-related genes"
)
dev.off()

# 2C. immune score boxplot
immune_score <- colMeans(mat_vsd[immune_idx, , drop = FALSE], na.rm = TRUE)
immune_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  immune_score = immune_score
)

p2C <- ggplot(immune_df, aes(x = group, y = immune_score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Figure 2C. Immune signature score") +
  xlab("") +
  ylab("Mean VST expression")

ggsave("Figures/Figure2C_immune_score_boxplot.png", p2C, width = 5, height = 4.5, dpi = 300)

# 2D. responder vs other volcano with immune labels
volcano_ro <- as.data.frame(res_RO)
volcano_ro$gene_id <- rownames(volcano_ro)
volcano_ro$gene_symbol <- sub(".*_", "", volcano_ro$gene_id)
volcano_ro <- volcano_ro[!is.na(volcano_ro$padj) & !is.na(volcano_ro$log2FoldChange), ]
volcano_ro$neglog10padj <- -log10(volcano_ro$padj)

volcano_ro$category <- "NS"
volcano_ro$category[volcano_ro$padj < 0.05 & volcano_ro$log2FoldChange > 1] <- "Up"
volcano_ro$category[volcano_ro$padj < 0.05 & volcano_ro$log2FoldChange < -1] <- "Down"

immune_label_df <- volcano_ro[volcano_ro$gene_symbol %in% immune_genes, ]

p2D <- ggplot(volcano_ro, aes(x = log2FoldChange, y = neglog10padj, color = category)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(
    data = immune_label_df,
    aes(label = gene_symbol),
    size = 3,
    vjust = -0.4,
    show.legend = FALSE
  ) +
  theme_bw() +
  ggtitle("Figure 2D. Volcano plot (Responder vs Other, immune genes)") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")

ggsave("Figures/Figure2D_volcano_immune_RO.png", p2D, width = 7, height = 5.5, dpi = 300)

# =========================================================
# Figure 3
# Melanogenesis program
# =========================================================

# 3A. melanogenesis heatmap
mel_idx <- which(gene_symbols_all %in% melanogenesis_genes)
mat_mel <- mat_vsd[mel_idx, , drop = FALSE]
rownames(mat_mel) <- make.unique(sub(".*_", "", rownames(mat_mel)))
mat_mel_scaled <- t(scale(t(mat_mel)))
mat_mel_scaled <- mat_mel_scaled[complete.cases(mat_mel_scaled), , drop = FALSE]

png("Figures/Figure3A_melanogenesis_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  mat_mel_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Figure 3A. Melanogenesis-related genes"
)
dev.off()

# 3B. melanogenesis score boxplot
mel_score <- colMeans(mat_vsd[mel_idx, , drop = FALSE], na.rm = TRUE)
mel_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  melanogenesis_score = mel_score
)

p3B <- ggplot(mel_df, aes(x = group, y = melanogenesis_score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Figure 3B. Melanogenesis signature score") +
  xlab("") +
  ylab("Mean VST expression")

ggsave("Figures/Figure3B_melanogenesis_score_boxplot.png", p3B, width = 5, height = 4.5, dpi = 300)

# 3C. responder vs other volcano with melanogenesis labels
mel_label_df <- volcano_ro[volcano_ro$gene_symbol %in% melanogenesis_genes, ]

p3C <- ggplot(volcano_ro, aes(x = log2FoldChange, y = neglog10padj, color = category)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(
    data = mel_label_df,
    aes(label = gene_symbol),
    size = 3,
    vjust = -0.4,
    show.legend = FALSE
  ) +
  theme_bw() +
  ggtitle("Figure 3C. Volcano plot (Responder vs Other, melanogenesis genes)") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")

ggsave("Figures/Figure3C_volcano_melanogenesis_RO.png", p3C, width = 7, height = 5.5, dpi = 300)

# 3D. mean score summary plot
mel_mean_df <- aggregate(melanogenesis_score ~ group, mel_df, mean)

p3D <- ggplot(mel_mean_df, aes(x = group, y = melanogenesis_score, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  ggtitle("Figure 3D. Mean melanogenesis score") +
  xlab("") +
  ylab("Mean score")

ggsave("Figures/Figure3D_mean_melanogenesis_score.png", p3D, width = 5, height = 4.5, dpi = 300)

# =========================================================
# Figure 4
# CTLA4 vs melanogenesis
# =========================================================

# 4A/B. scatter + regression
ctla4_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  Ctla4 = mat_vsd[ctla4_idx[1], ],
  melanogenesis_score = mel_score
)

p4A <- ggplot(ctla4_df, aes(x = Ctla4, y = melanogenesis_score, color = group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Figure 4A. Ctla4 vs melanogenesis score") +
  xlab("Ctla4 expression (VST)") +
  ylab("Melanogenesis score")

ggsave("Figures/Figure4A_Ctla4_vs_melanogenesis.png", p4A, width = 6, height = 5, dpi = 300)

# 4C. side-by-side boxplots
p4C1 <- p2A + ggtitle("Ctla4 expression")
p4C2 <- p3B + ggtitle("Melanogenesis score")

png("Figures/Figure4C_Ctla4_and_melanogenesis_boxplots.png", width = 1600, height = 700, res = 150)
grid.arrange(p4C1, p4C2, ncol = 2)
dev.off()

# 4D. simple conceptual summary plot
concept_df <- data.frame(
  variable = c("Ctla4", "Melanogenesis"),
  Control = c(mean(df_ctla4$expression[df_ctla4$group == "Control"]),
              mean(mel_df$melanogenesis_score[mel_df$group == "Control"])),
  Other = c(mean(df_ctla4$expression[df_ctla4$group == "Other"]),
            mean(mel_df$melanogenesis_score[mel_df$group == "Other"])),
  Responder = c(mean(df_ctla4$expression[df_ctla4$group == "Responder"]),
                mean(mel_df$melanogenesis_score[mel_df$group == "Responder"]))
)

concept_long <- melt(concept_df, id.vars = "variable", variable.name = "group", value.name = "value")

p4D <- ggplot(concept_long, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  ggtitle("Figure 4D. Opposite trends of Ctla4 and melanogenesis") +
  xlab("") +
  ylab("Mean value")

ggsave("Figures/Figure4D_concept_barplot.png", p4D, width = 6.5, height = 5, dpi = 300)

# =========================================================
# Figure 5
# Mechanistic summary / integrated state model
# =========================================================

# 5A. top variable genes heatmap (already biologically informative)
top_var <- apply(mat_vsd, 1, var)
top_genes <- names(sort(top_var, decreasing = TRUE))[1:50]
mat_top <- mat_vsd[top_genes, , drop = FALSE]
rownames(mat_top) <- make.unique(sub(".*_", "", rownames(mat_top)))
mat_top_scaled <- t(scale(t(mat_top)))
mat_top_scaled <- mat_top_scaled[complete.cases(mat_top_scaled), , drop = FALSE]

png("Figures/Figure5A_top50_variable_genes_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  mat_top_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  main = "Figure 5A. Top 50 variable genes"
)
dev.off()

# 5B. responder vs other selected genes barplot
focus_ro <- res_RO_df[res_RO_df$gene_symbol %in% c("Ctla4", "Dct", "Tyr", "Pmel", "Mlana", "Stat1", "Cxcl9"), ]
focus_ro_plot <- focus_ro[, c("gene_symbol", "log2FoldChange")]
focus_ro_plot$gene_symbol <- factor(focus_ro_plot$gene_symbol, levels = focus_ro_plot$gene_symbol)

p5B <- ggplot(focus_ro_plot, aes(x = gene_symbol, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  ggtitle("Figure 5B. Key genes in Responder vs Other") +
  xlab("") +
  ylab("log2 fold change")

ggsave("Figures/Figure5B_key_genes_RO.png", p5B, width = 6.5, height = 5, dpi = 300)

# 5C. immune vs melanogenesis score scatter
score_df <- data.frame(
  sample = colnames(vsd),
  group = meta$group,
  immune_score = immune_score,
  melanogenesis_score = mel_score
)

p5C <- ggplot(score_df, aes(x = immune_score, y = melanogenesis_score, color = group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  ggtitle("Figure 5C. Immune vs melanogenesis score") +
  xlab("Immune score") +
  ylab("Melanogenesis score")

ggsave("Figures/Figure5C_immune_vs_melanogenesis.png", p5C, width = 6, height = 5, dpi = 300)

# 5D. summary table as CSV for manuscript
summary_stats <- list(
  DEG_counts = deg_counts,
  Melanogenesis_group_mean = mel_mean_df,
  Ctla4_melanogenesis_correlation = data.frame(
    metric = "Pearson_correlation",
    value = cor(ctla4_df$Ctla4, ctla4_df$melanogenesis_score, method = "pearson")
  )
)

write.csv(deg_counts, "Figures/Figure5D_DEG_counts_summary.csv", row.names = FALSE)
write.csv(mel_mean_df, "Figures/Figure5D_melanogenesis_group_mean.csv", row.names = FALSE)
write.csv(summary_stats$Ctla4_melanogenesis_correlation,
          "Figures/Figure5D_Ctla4_melanogenesis_correlation.csv",
          row.names = FALSE)

cat("All figure files saved in ./Figures\n")