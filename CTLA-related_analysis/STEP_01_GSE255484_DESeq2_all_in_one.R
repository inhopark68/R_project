# =========================================================
# GSE255484 DESeq2 통합 분석 스크립트
# - raw counts 불러오기
# - meta 매칭
# - group 생성
# - DESeq2 실행
# - gene symbol annotation
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
library(AnnotationDbi)
library(pheatmap)

# 종 선택
# mouse / human 중 하나
species_db <- "mouse"

if (species_db == "mouse") {
  library(org.Mm.eg.db)
  anno_db <- org.Mm.eg.db
} else if (species_db == "human") {
  library(org.Hs.eg.db)
  anno_db <- org.Hs.eg.db
} else {
  stop("species_db must be either 'mouse' or 'human'")
}

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
# stable + non-responder -> Non-responder
#
# 중요:
# "non-responder" 안에 "responder" 문자열이 포함되므로
# 반드시 non-responder/stable 을 responder보다 먼저 판별해야 함
# -----------------------------
title_clean <- trimws(meta$title)

meta$group <- ifelse(
  grepl("\\bcontrol\\b", title_clean, ignore.case = TRUE),
  "Control",
  ifelse(
    grepl("non[- ]?responder|\\bstable\\b", title_clean, ignore.case = TRUE),
    "Non-responder",
    ifelse(
      grepl("\\bresponder\\b", title_clean, ignore.case = TRUE),
      "Responder",
      NA_character_
    )
  )
)

meta$group <- factor(
  meta$group,
  levels = c("Control", "Non-responder", "Responder")
)

cat("group table:\n")
print(table(meta$group, useNA = "ifany"))

if (any(is.na(meta$group))) {
  warning("Some samples could not be assigned to a group. Check meta$title values.")
}

# DESeq2에는 NA group이 들어가면 안 되므로 제거
keep_idx <- !is.na(meta$group)
meta <- meta[keep_idx, , drop = FALSE]
expr <- expr[, keep_idx, drop = FALSE]

cat("samples retained after group assignment:", ncol(expr), "\n")

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
res_NC <- results(dds, contrast = c("group", "Non-responder", "Control"))
res_RN <- results(dds, contrast = c("group", "Responder", "Non-responder"))

# -----------------------------
# 7) gene annotation 함수
# -----------------------------
add_gene_annotation <- function(res_obj, anno_db) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_id_clean <- sub("\\..*$", "", res_df$gene_id)

  res_df$gene_symbol <- mapIds(
    anno_db,
    keys = res_df$gene_id_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

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

  res_df$gene_label <- ifelse(
    is.na(res_df$gene_symbol) | res_df$gene_symbol == "",
    res_df$gene_id,
    res_df$gene_symbol
  )

  res_df <- res_df[order(res_df$padj), ]
  res_df
}

# annotation 추가
res_RC_df <- add_gene_annotation(res_RC, anno_db)
res_NC_df <- add_gene_annotation(res_NC, anno_db)
res_RN_df <- add_gene_annotation(res_RN, anno_db)

# -----------------------------
# 8) 유의 유전자 추출
# -----------------------------
sig_RC <- subset(res_RC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_NC <- subset(res_NC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_RN <- subset(res_RN_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# -----------------------------
# 9) 결과 저장
# -----------------------------
write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(res_NC_df, "DEG_NonResponder_vs_Control.csv", row.names = FALSE)
write.csv(res_RN_df, "DEG_Responder_vs_NonResponder.csv", row.names = FALSE)

write.csv(sig_RC, "SIG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(sig_NC, "SIG_NonResponder_vs_Control.csv", row.names = FALSE)
write.csv(sig_RN, "SIG_Responder_vs_NonResponder.csv", row.names = FALSE)

# -----------------------------
# 10) MA plot 저장
# -----------------------------
png("MAplot_Responder_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(res_RC, main = "MA plot: Responder vs Control", ylim = c(-5, 5))
dev.off()

png("MAplot_NonResponder_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(res_NC, main = "MA plot: Non-responder vs Control", ylim = c(-5, 5))
dev.off()

png("MAplot_Responder_vs_NonResponder.png", width = 1200, height = 1000, res = 150)
plotMA(res_RN, main = "MA plot: Responder vs Non-responder", ylim = c(-5, 5))
dev.off()

# -----------------------------
# 11) VST 변환
# -----------------------------
vsd <- vst(dds, blind = FALSE)

# -----------------------------
# 12) Sample distance heatmap 저장
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
# 13) normalized counts 저장
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)

norm_counts_df$gene_id <- rownames(norm_counts_df)
norm_counts_df$gene_id_clean <- sub("\\..*$", "", norm_counts_df$gene_id)

norm_counts_df$gene_symbol <- mapIds(
  anno_db,
  keys = norm_counts_df$gene_id_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

norm_counts_df$gene_name <- mapIds(
  anno_db,
  keys = norm_counts_df$gene_id_clean,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

norm_counts_df <- norm_counts_df[, c(
  "gene_id",
  "gene_id_clean",
  "gene_symbol",
  "gene_name",
  setdiff(colnames(norm_counts_df), c("gene_id", "gene_id_clean", "gene_symbol", "gene_name"))
)]

write.csv(norm_counts_df, "Normalized_counts.csv", row.names = FALSE)

# -----------------------------
# 14) top variable genes heatmap 저장
# -----------------------------
mat <- assay(vsd)
gene_var <- apply(mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, nrow(mat))]
mat_top <- mat[top_genes, , drop = FALSE]

top_genes_clean <- sub("\\..*$", "", rownames(mat_top))
top_symbols <- mapIds(
  anno_db,
  keys = top_genes_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

row_labels <- ifelse(is.na(top_symbols) | top_symbols == "", rownames(mat_top), top_symbols)
rownames(mat_top) <- make.unique(row_labels)

mat_top_scaled <- t(scale(t(mat_top)))

annotation_col <- data.frame(group = colData(dds)$group)
rownames(annotation_col) <- colnames(dds)

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
# 15) 완료 메시지
# -----------------------------
cat("\n===== 분석 완료 =====\n")
cat("유의 유전자 수\n")
cat("Responder vs Control:", nrow(sig_RC), "\n")
cat("Non-responder vs Control:", nrow(sig_NC), "\n")
cat("Responder vs Non-responder:", nrow(sig_RN), "\n")

cat("\n저장된 파일\n")
cat("- DEG_Responder_vs_Control.csv\n")
cat("- DEG_NonResponder_vs_Control.csv\n")
cat("- DEG_Responder_vs_NonResponder.csv\n")
cat("- SIG_Responder_vs_Control.csv\n")
cat("- SIG_NonResponder_vs_Control.csv\n")
cat("- SIG_Responder_vs_NonResponder.csv\n")
cat("- MAplot_Responder_vs_Control.png\n")
cat("- MAplot_NonResponder_vs_Control.png\n")
cat("- MAplot_Responder_vs_NonResponder.png\n")
cat("- Sample_distance_heatmap.png\n")
cat("- Normalized_counts.csv\n")
cat("- Heatmap_top50_variable_genes.png\n")