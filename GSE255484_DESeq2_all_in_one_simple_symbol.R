# =========================================================
# GSE255484 DESeq2 통합 분석 스크립트
# - raw counts 불러오기
# - meta 매칭
# - group 생성
# - DESeq2 실행
# - gene symbol 추출 (언더바 뒤 문자열)
# - DEG 저장
# - MA plot 저장
# - sample distance heatmap 저장
# - normalized counts 저장
# =========================================================

# -----------------------------
# 0) 패키지 로드
# -----------------------------
library(GEOquery)
library(DESeq2)
library(pheatmap)

# -----------------------------
# 1) GEO 메타데이터 불러오기
# -----------------------------
gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

cat("meta dim:", dim(meta), "\n")
cat("meta preview:\n")
print(head(meta$title, 10))

# -----------------------------
# 2) raw count 파일 다운로드 및 읽기
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
# 3) 샘플 매칭
# expr 열 이름 예: 4_166__EPG
# meta$title 예: tumor, anti-CTLA-4, responder, id 166
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
# 4) group 생성
# responder -> Responder
# control -> Control
# stable + non-responder -> Other
# -----------------------------
meta$group <- ifelse(
  grepl("responder", meta$title, ignore.case = TRUE),
  "Responder",
  ifelse(
    grepl("control", meta$title, ignore.case = TRUE),
    "Control",
    "Other"
  )
)

meta$group <- factor(meta$group, levels = c("Control", "Other", "Responder"))

cat("group table:\n")
print(table(meta$group))

# -----------------------------
# 5) DESeq2 객체 생성 및 실행
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
# 6) DEG 결과 추출
# -----------------------------
res_RC <- results(dds, contrast = c("group", "Responder", "Control"))
res_OC <- results(dds, contrast = c("group", "Other", "Control"))
res_RO <- results(dds, contrast = c("group", "Responder", "Other"))

# -----------------------------
# 7) gene symbol 추출 함수
# 예: ENSMUSG00000000001.4_Gnai3 -> Gnai3
# -----------------------------
add_gene_symbol <- function(res_obj) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_symbol <- sub(".*_", "", res_df$gene_id)
  res_df <- res_df[order(res_df$padj), ]
  res_df
}

res_RC_df <- add_gene_symbol(res_RC)
res_OC_df <- add_gene_symbol(res_OC)
res_RO_df <- add_gene_symbol(res_RO)

# -----------------------------
# 8) 유의 유전자 추출
# -----------------------------
sig_RC <- subset(res_RC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_OC <- subset(res_OC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_RO <- subset(res_RO_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# -----------------------------
# 9) 컬럼 순서 정리
# -----------------------------
move_gene_cols_front <- function(df) {
  df[, c("gene_id", "gene_symbol", setdiff(colnames(df), c("gene_id", "gene_symbol")))]
}

res_RC_df <- move_gene_cols_front(res_RC_df)
res_OC_df <- move_gene_cols_front(res_OC_df)
res_RO_df <- move_gene_cols_front(res_RO_df)

sig_RC <- move_gene_cols_front(sig_RC)
sig_OC <- move_gene_cols_front(sig_OC)
sig_RO <- move_gene_cols_front(sig_RO)

# -----------------------------
# 10) 결과 저장
# -----------------------------
write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(res_OC_df, "DEG_Other_vs_Control.csv", row.names = FALSE)
write.csv(res_RO_df, "DEG_Responder_vs_Other.csv", row.names = FALSE)

write.csv(sig_RC, "SIG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(sig_OC, "SIG_Other_vs_Control.csv", row.names = FALSE)
write.csv(sig_RO, "SIG_Responder_vs_Other.csv", row.names = FALSE)

# -----------------------------
# 11) MA plot 저장
# -----------------------------
png("MAplot_Responder_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(res_RC, main = "MA plot: Responder vs Control", ylim = c(-5, 5))
dev.off()

png("MAplot_Other_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(res_OC, main = "MA plot: Other vs Control", ylim = c(-5, 5))
dev.off()

png("MAplot_Responder_vs_Other.png", width = 1200, height = 1000, res = 150)
plotMA(res_RO, main = "MA plot: Responder vs Other", ylim = c(-5, 5))
dev.off()

# -----------------------------
# 12) VST 변환
# -----------------------------
vsd <- vst(dds, blind = FALSE)

# -----------------------------
# 13) Sample distance heatmap 저장
# -----------------------------
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

sample_labels <- paste(colnames(vsd), colData(vsd)$group, sep = "_")
rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels

png("Sample_distance_heatmap.png", width = 1400, height = 1200, res = 150)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample-to-sample distances"
)
dev.off()

# -----------------------------
# 14) normalized counts 저장
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)

norm_counts_df$gene_id <- rownames(norm_counts_df)
norm_counts_df$gene_symbol <- sub(".*_", "", norm_counts_df$gene_id)

norm_counts_df <- norm_counts_df[, c(
  "gene_id",
  "gene_symbol",
  setdiff(colnames(norm_counts_df), c("gene_id", "gene_symbol"))
)]

write.csv(norm_counts_df, "Normalized_counts.csv", row.names = FALSE)

# -----------------------------
# 15) top variable genes heatmap 저장
# -----------------------------
mat <- assay(vsd)
gene_var <- apply(mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, nrow(mat))]
mat_top <- mat[top_genes, , drop = FALSE]

row_labels <- sub(".*_", "", rownames(mat_top))
rownames(mat_top) <- make.unique(row_labels)

mat_top_scaled <- t(scale(t(mat_top)))

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

png("Heatmap_top50_variable_genes.png", width = 1200, height = 1000, res = 150)
pheatmap(
  mat_top_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  main = "Top 50 variable genes"
)
dev.off()

# -----------------------------
# 16) 완료 메시지
# -----------------------------
cat("\n===== 분석 완료 =====\n")
cat("유의 유전자 수\n")
cat("Responder vs Control:", nrow(sig_RC), "\n")
cat("Other vs Control:", nrow(sig_OC), "\n")
cat("Responder vs Other:", nrow(sig_RO), "\n")

cat("\n저장된 파일\n")
cat("- DEG_Responder_vs_Control.csv\n")
cat("- DEG_Other_vs_Control.csv\n")
cat("- DEG_Responder_vs_Other.csv\n")
cat("- SIG_Responder_vs_Control.csv\n")
cat("- SIG_Other_vs_Control.csv\n")
cat("- SIG_Responder_vs_Other.csv\n")
cat("- MAplot_Responder_vs_Control.png\n")
cat("- MAplot_Other_vs_Control.png\n")
cat("- MAplot_Responder_vs_Other.png\n")
cat("- Sample_distance_heatmap.png\n")
cat("- Normalized_counts.csv\n")
cat("- Heatmap_top50_variable_genes.png\n")