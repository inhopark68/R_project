# =========================================================
# GSE255484 integrated analysis + Figure 1~5 generation
# Final integrated version
# - group split: Control / Stable / Non-responder / Responder
# - robust group parsing to avoid NA labels
# - DESeq2 analysis
# - volcano plots with immune / melanogenesis gene labels
# - Figure 1~5 generation
# =========================================================

# -----------------------------
# 0) package load
# -----------------------------
library(GEOquery)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(dplyr)
library(tidyr)

# -----------------------------
# 1) output directories
# -----------------------------
dir.create("Results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

# -----------------------------
# 2) annotation DB
# -----------------------------
anno_db <- org.Mm.eg.db

# -----------------------------
# 3) GEO metadata load
# -----------------------------
gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

cat("meta dim:", dim(meta), "\n")
cat("meta preview:\n")
print(head(meta$title, 10))

# -----------------------------
# 4) raw count load
# -----------------------------
getGEOSuppFiles("GSE255484", makeDirectory = TRUE)

count_file <- "GSE255484/GSE255484_raw_counts.txt.gz"

count_df <- read.delim(
  count_file,
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(count_df)
mode(expr) <- "numeric"

cat("expr dim:", dim(expr), "\n")
cat("expr colnames preview:\n")
print(head(colnames(expr), 10))

# -----------------------------
# 5) sample matching
# expr colnames example: 4_166__EPG
# meta$title example: tumor, anti-CTLA-4, responder, id 166
# -----------------------------
expr_id <- sub("^[0-9]+_([0-9]+)__EPG$", "\\1", colnames(expr))
meta$id_in_title <- sub(".*id[[:space:]]+([0-9]+).*", "\\1", meta$title)

meta <- meta[match(expr_id, meta$id_in_title), ]
rownames(meta) <- colnames(expr)

cat("sample matched:", all(colnames(expr) == rownames(meta)), "\n")
cat("NA in matched meta:", sum(is.na(meta$id_in_title)), "\n")

if (!all(colnames(expr) == rownames(meta))) {
  stop("Sample matching failed: colnames(expr) != rownames(meta)")
}

# -----------------------------
# 6) robust group assignment
# WHY 'NA' happened:
# samples not matching the old rules were assigned NA.
# Here we assign Unclassified first, then exclude them from DESeq2/Figures.
# -----------------------------
title_clean <- trimws(meta$title)

meta$group <- case_when(
  grepl("\\bcontrol\\b|\\buntreated\\b|\\bvehicle\\b|\\bigg\\b", title_clean, ignore.case = TRUE) ~ "Control",
  grepl("non[- ]?responder|non[- ]?response|\\bnr\\b", title_clean, ignore.case = TRUE) ~ "Non-responder",
  grepl("\\bstable\\b|stable disease|\\bsd\\b", title_clean, ignore.case = TRUE) ~ "Stable",
  grepl("\\bresponder\\b|\\bresponse\\b", title_clean, ignore.case = TRUE) ~ "Responder",
  TRUE ~ "Unclassified"
)

meta$group <- factor(
  meta$group,
  levels = c("Control", "Stable", "Non-responder", "Responder", "Unclassified")
)

cat("group table before filtering:\n")
print(table(meta$group, useNA = "ifany"))

if (any(meta$group == "Unclassified")) {
  cat("\nUnclassified sample titles:\n")
  print(unique(meta$title[meta$group == "Unclassified"]))
}

# Remove unclassified samples for downstream analysis/figures
keep_idx <- meta$group != "Unclassified" & !is.na(meta$group)
meta <- meta[keep_idx, , drop = FALSE]
expr <- expr[, keep_idx, drop = FALSE]

meta$group <- factor(
  as.character(meta$group),
  levels = c("Control", "Stable", "Non-responder", "Responder")
)

cat("group table after filtering:\n")
print(table(meta$group, useNA = "ifany"))
cat("samples retained after group assignment:", ncol(expr), "\n")

# -----------------------------
# 7) DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = meta,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

cat("dds dim after filtering:", dim(dds), "\n")

# -----------------------------
# 8) contrasts
# -----------------------------
res_RC <- results(dds, contrast = c("group", "Responder", "Control"))
res_SC <- results(dds, contrast = c("group", "Stable", "Control"))
res_NC <- results(dds, contrast = c("group", "Non-responder", "Control"))
res_RS <- results(dds, contrast = c("group", "Responder", "Stable"))
res_RN <- results(dds, contrast = c("group", "Responder", "Non-responder"))
res_NS <- results(dds, contrast = c("group", "Non-responder", "Stable"))

# -----------------------------
# 9) result annotation helpers
# gene format example: ENSMUSG...._Gnai3
# -----------------------------
add_gene_annotation <- function(res_obj, anno_db) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_symbol <- sub(".*_", "", res_df$gene_id)
  res_df$gene_id_clean <- sub("\\..*$", "", sub("_.*$", "", res_df$gene_id))

  suppressWarnings({
    res_df$entrez_id <- mapIds(
      anno_db,
      keys = res_df$gene_id_clean,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

    res_df$gene_name <- mapIds(
      anno_db,
      keys = res_df$gene_id_clean,
      column = "GENENAME",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  })

  res_df <- res_df[order(res_df$padj), ]
  res_df
}

move_gene_cols_front <- function(df) {
  df[, c(
    "gene_id", "gene_symbol",
    setdiff(colnames(df), c("gene_id", "gene_symbol"))
  )]
}

get_sig <- function(df) {
  subset(df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
}

res_RC_df <- move_gene_cols_front(add_gene_annotation(res_RC, anno_db))
res_SC_df <- move_gene_cols_front(add_gene_annotation(res_SC, anno_db))
res_NC_df <- move_gene_cols_front(add_gene_annotation(res_NC, anno_db))
res_RS_df <- move_gene_cols_front(add_gene_annotation(res_RS, anno_db))
res_RN_df <- move_gene_cols_front(add_gene_annotation(res_RN, anno_db))
res_NS_df <- move_gene_cols_front(add_gene_annotation(res_NS, anno_db))

sig_RC <- move_gene_cols_front(get_sig(res_RC_df))
sig_SC <- move_gene_cols_front(get_sig(res_SC_df))
sig_NC <- move_gene_cols_front(get_sig(res_NC_df))
sig_RS <- move_gene_cols_front(get_sig(res_RS_df))
sig_RN <- move_gene_cols_front(get_sig(res_RN_df))
sig_NS <- move_gene_cols_front(get_sig(res_NS_df))

# -----------------------------
# 10) save DEG results
# -----------------------------
write.csv(res_RC_df, "Results/DEG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(res_SC_df, "Results/DEG_Stable_vs_Control.csv", row.names = FALSE)
write.csv(res_NC_df, "Results/DEG_NonResponder_vs_Control.csv", row.names = FALSE)
write.csv(res_RS_df, "Results/DEG_Responder_vs_Stable.csv", row.names = FALSE)
write.csv(res_RN_df, "Results/DEG_Responder_vs_NonResponder.csv", row.names = FALSE)
write.csv(res_NS_df, "Results/DEG_NonResponder_vs_Stable.csv", row.names = FALSE)

write.csv(sig_RC, "Results/SIG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(sig_SC, "Results/SIG_Stable_vs_Control.csv", row.names = FALSE)
write.csv(sig_NC, "Results/SIG_NonResponder_vs_Control.csv", row.names = FALSE)
write.csv(sig_RS, "Results/SIG_Responder_vs_Stable.csv", row.names = FALSE)
write.csv(sig_RN, "Results/SIG_Responder_vs_NonResponder.csv", row.names = FALSE)
write.csv(sig_NS, "Results/SIG_NonResponder_vs_Stable.csv", row.names = FALSE)

# -----------------------------
# 11) MA plots
# -----------------------------
save_ma_plot <- function(res_obj, filename, title_text) {
  png(filename, width = 1200, height = 1000, res = 150)
  plotMA(res_obj, main = title_text, ylim = c(-5, 5))
  dev.off()
}

save_ma_plot(res_RC, "Results/MAplot_Responder_vs_Control.png", "MA plot: Responder vs Control")
save_ma_plot(res_SC, "Results/MAplot_Stable_vs_Control.png", "MA plot: Stable vs Control")
save_ma_plot(res_NC, "Results/MAplot_NonResponder_vs_Control.png", "MA plot: Non-responder vs Control")
save_ma_plot(res_RS, "Results/MAplot_Responder_vs_Stable.png", "MA plot: Responder vs Stable")
save_ma_plot(res_RN, "Results/MAplot_Responder_vs_NonResponder.png", "MA plot: Responder vs Non-responder")
save_ma_plot(res_NS, "Results/MAplot_NonResponder_vs_Stable.png", "MA plot: Non-responder vs Stable")

# -----------------------------
# 12) VST
# -----------------------------
vsd <- vst(dds, blind = FALSE)
mat_vsd <- assay(vsd)
plot_group <- factor(
  colData(vsd)$group,
  levels = c("Control", "Stable", "Non-responder", "Responder")
)

# -----------------------------
# 13) save normalized counts
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene_id <- rownames(norm_counts_df)
norm_counts_df$gene_symbol <- sub(".*_", "", norm_counts_df$gene_id)
norm_counts_df <- norm_counts_df[, c(
  "gene_id", "gene_symbol",
  setdiff(colnames(norm_counts_df), c("gene_id", "gene_symbol"))
)]

write.csv(norm_counts_df, "Results/Normalized_counts.csv", row.names = FALSE)

# -----------------------------
# 14) gene sets
# -----------------------------
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
annotation_col <- data.frame(group = plot_group)
rownames(annotation_col) <- colnames(vsd)

# -----------------------------
# 15) helper functions for figure code
# -----------------------------
make_volcano_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  df$gene_id <- rownames(df)
  df$gene_symbol <- sub(".*_", "", df$gene_id)
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), , drop = FALSE]
  df$neglog10padj <- -log10(df$padj)

  df$category <- "NS"
  df$category[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
  df$category[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"
  df
}

make_sig_label_df <- function(volcano_df, genes) {
  volcano_df[
    volcano_df$gene_symbol %in% genes &
      volcano_df$padj < 0.05 &
      abs(volcano_df$log2FoldChange) > 1,
    ,
    drop = FALSE
  ]
}

plot_volcano_with_labels <- function(volcano_df, label_df, title_text) {
  ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10padj, color = category)) +
    geom_point(alpha = 0.6, size = 1.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
      data = label_df,
      aes(label = gene_symbol),
      size = 3.5,
      color = "black",
      box.padding = 0.4,
      point.padding = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    theme_bw() +
    ggtitle(title_text) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value")
}

# Volcano dfs
volcano_RN <- make_volcano_df(res_RN)
volcano_RS <- make_volcano_df(res_RS)
volcano_NS <- make_volcano_df(res_NS)

# -----------------------------
# Figure 1
# -----------------------------

# Figure 1A
df_group <- as.data.frame(table(plot_group))
colnames(df_group) <- c("group", "count")

p1A <- ggplot(df_group, aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  ggtitle("Figure 1A. Sample groups") +
  xlab("") +
  ylab("Number of samples")

ggsave("Figures/Figure1A_group_counts.png", p1A, width = 6, height = 4, dpi = 300)

# Figure 1B
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p1B <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("Figure 1B. PCA")

ggsave("Figures/Figure1B_PCA.png", p1B, width = 6, height = 5, dpi = 300)

# Figure 1C
sample_dists <- dist(t(mat_vsd))
sample_dist_matrix <- as.matrix(sample_dists)
sample_labels <- paste(colnames(vsd), plot_group, sep = "_")
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

# Figure 1D
deg_counts <- data.frame(
  comparison = c(
    "Responder vs Control",
    "Stable vs Control",
    "Non-responder vs Control",
    "Responder vs Stable",
    "Responder vs Non-responder",
    "Non-responder vs Stable"
  ),
  n_sig = c(
    nrow(sig_RC),
    nrow(sig_SC),
    nrow(sig_NC),
    nrow(sig_RS),
    nrow(sig_RN),
    nrow(sig_NS)
  )
)

p1D <- ggplot(deg_counts, aes(x = comparison, y = n_sig, fill = comparison)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  ggtitle("Figure 1D. Number of significant DEGs") +
  xlab("") +
  ylab("DEG count") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("Figures/Figure1D_DEG_counts.png", p1D, width = 8, height = 5, dpi = 300)

# -----------------------------
# Figure 2
# -----------------------------

# Figure 2A
ctla4_idx <- which(gene_symbols_all == "Ctla4")
df_ctla4 <- data.frame(
  expression = mat_vsd[ctla4_idx[1], ],
  group = plot_group
)

p2A <- ggplot(df_ctla4, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Figure 2A. Ctla4 expression") +
  xlab("") +
  ylab("VST expression")

ggsave("Figures/Figure2A_Ctla4_boxplot.png", p2A, width = 6, height = 4.5, dpi = 300)

# Figure 2B
immune_idx <- which(gene_symbols_all %in% immune_genes)
mat_immune <- mat_vsd[immune_idx, , drop = FALSE]
rownames(mat_immune) <- make.unique(sub(".*_", "", rownames(mat_immune)))
mat_immune_scaled <- t(scale(t(mat_immune)))
mat_immune_scaled <- mat_immune_scaled[complete.cases(mat_immune_scaled), , drop = FALSE]

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

# Figure 2C
immune_score <- colMeans(mat_vsd[immune_idx, , drop = FALSE], na.rm = TRUE)
immune_df <- data.frame(
  sample = colnames(vsd),
  group = plot_group,
  immune_score = immune_score
)

p2C <- ggplot(immune_df, aes(x = group, y = immune_score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Figure 2C. Immune signature score") +
  xlab("") +
  ylab("Mean VST expression")

ggsave("Figures/Figure2C_immune_score_boxplot.png", p2C, width = 6, height = 4.5, dpi = 300)

# Figure 2D
immune_label_df_RN <- make_sig_label_df(volcano_RN, immune_genes)

p2D <- plot_volcano_with_labels(
  volcano_df = volcano_RN,
  label_df = immune_label_df_RN,
  title_text = "Figure 2D. Volcano plot (Responder vs Non-responder, immune genes)"
)

ggsave("Figures/Figure2D_volcano_immune_RN.png", p2D, width = 7.5, height = 5.5, dpi = 300)

# -----------------------------
# Figure 3
# -----------------------------

# Figure 3A
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

# Figure 3B
mel_score <- colMeans(mat_vsd[mel_idx, , drop = FALSE], na.rm = TRUE)
mel_df <- data.frame(
  sample = colnames(vsd),
  group = plot_group,
  melanogenesis_score = mel_score
)

p3B <- ggplot(mel_df, aes(x = group, y = melanogenesis_score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  ggtitle("Figure 3B. Melanogenesis signature score") +
  xlab("") +
  ylab("Mean VST expression")

ggsave("Figures/Figure3B_melanogenesis_score_boxplot.png", p3B, width = 6, height = 4.5, dpi = 300)

# Figure 3C
mel_label_df_RN <- make_sig_label_df(volcano_RN, melanogenesis_genes)

p3C <- plot_volcano_with_labels(
  volcano_df = volcano_RN,
  label_df = mel_label_df_RN,
  title_text = "Figure 3C. Volcano plot (Responder vs Non-responder, melanogenesis genes)"
)

ggsave("Figures/Figure3C_volcano_melanogenesis_RN.png", p3C, width = 7.5, height = 5.5, dpi = 300)

# Figure 3D
mel_mean_df <- aggregate(melanogenesis_score ~ group, mel_df, mean)

p3D <- ggplot(mel_mean_df, aes(x = group, y = melanogenesis_score, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  ggtitle("Figure 3D. Mean melanogenesis score") +
  xlab("") +
  ylab("Mean score")

ggsave("Figures/Figure3D_mean_melanogenesis_score.png", p3D, width = 6, height = 4.5, dpi = 300)

# -----------------------------
# Figure 4
# -----------------------------

# Figure 4A
ctla4_df <- data.frame(
  sample = colnames(vsd),
  group = plot_group,
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

# Figure 4B
p4B1 <- p2A + ggtitle("Ctla4 expression")
p4B2 <- p3B + ggtitle("Melanogenesis score")

png("Figures/Figure4B_Ctla4_and_melanogenesis_boxplots.png", width = 1600, height = 700, res = 150)
grid.arrange(p4B1, p4B2, ncol = 2)
dev.off()

# Figure 4C
concept_df <- data.frame(
  variable = c("Ctla4", "Melanogenesis"),
  Control = c(
    mean(df_ctla4$expression[df_ctla4$group == "Control"]),
    mean(mel_df$melanogenesis_score[mel_df$group == "Control"])
  ),
  Stable = c(
    mean(df_ctla4$expression[df_ctla4$group == "Stable"]),
    mean(mel_df$melanogenesis_score[mel_df$group == "Stable"])
  ),
  `Non-responder` = c(
    mean(df_ctla4$expression[df_ctla4$group == "Non-responder"]),
    mean(mel_df$melanogenesis_score[mel_df$group == "Non-responder"])
  ),
  Responder = c(
    mean(df_ctla4$expression[df_ctla4$group == "Responder"]),
    mean(mel_df$melanogenesis_score[mel_df$group == "Responder"])
  ),
  check.names = FALSE
)

concept_long <- melt(concept_df, id.vars = "variable", variable.name = "group", value.name = "value")
concept_long$group <- factor(concept_long$group, levels = c("Control", "Stable", "Non-responder", "Responder"))

p4C <- ggplot(concept_long, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  ggtitle("Figure 4C. Opposite trends of Ctla4 and melanogenesis") +
  xlab("") +
  ylab("Mean value")

ggsave("Figures/Figure4C_concept_barplot.png", p4C, width = 7, height = 5, dpi = 300)

# Figure 4D
ctla4_mel_cor <- data.frame(
  metric = "Pearson_correlation",
  value = cor(ctla4_df$Ctla4, ctla4_df$melanogenesis_score, method = "pearson")
)

write.csv(ctla4_mel_cor, "Figures/Figure4D_Ctla4_melanogenesis_correlation.csv", row.names = FALSE)

# -----------------------------
# Figure 5
# -----------------------------

# Figure 5A
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

# Figure 5B
focus_rn <- res_RN_df[
  res_RN_df$gene_symbol %in% c("Ctla4", "Dct", "Tyr", "Pmel", "Mlana", "Stat1", "Cxcl9"),
]
focus_rn_plot <- focus_rn[, c("gene_symbol", "log2FoldChange")]
focus_rn_plot$gene_symbol <- factor(focus_rn_plot$gene_symbol, levels = focus_rn_plot$gene_symbol)

p5B <- ggplot(focus_rn_plot, aes(x = gene_symbol, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  ggtitle("Figure 5B. Key genes in Responder vs Non-responder") +
  xlab("") +
  ylab("log2 fold change")

ggsave("Figures/Figure5B_key_genes_RN.png", p5B, width = 6.5, height = 5, dpi = 300)

# Figure 5C
score_df <- data.frame(
  sample = colnames(vsd),
  group = plot_group,
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

# Figure 5D
write.csv(deg_counts, "Figures/Figure5D_DEG_counts_summary.csv", row.names = FALSE)
write.csv(mel_mean_df, "Figures/Figure5D_melanogenesis_group_mean.csv", row.names = FALSE)
write.csv(ctla4_mel_cor, "Figures/Figure5D_Ctla4_melanogenesis_correlation.csv", row.names = FALSE)

# -----------------------------
# 16) final summary
# -----------------------------
cat("\n===== Analysis complete =====\n")
cat("Significant DEG counts\n")
cat("Responder vs Control:", nrow(sig_RC), "\n")
cat("Stable vs Control:", nrow(sig_SC), "\n")
cat("Non-responder vs Control:", nrow(sig_NC), "\n")
cat("Responder vs Stable:", nrow(sig_RS), "\n")
cat("Responder vs Non-responder:", nrow(sig_RN), "\n")
cat("Non-responder vs Stable:", nrow(sig_NS), "\n")

cat("\nSaved result files in ./Results\n")
cat("Saved figure files in ./Figures\n")