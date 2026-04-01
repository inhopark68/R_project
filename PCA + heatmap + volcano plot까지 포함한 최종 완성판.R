# ============================================
# GSE255484 DESeq2 전체 분석 + PCA + Heatmap + Volcano
# ============================================

# 0) 패키지 로드
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# ggrepel은 있으면 라벨이 더 예쁘게 나옵니다.
# 없으면 자동으로 라벨 없이 진행합니다.
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# 1) GEO 메타데이터 불러오기
gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

# 2) raw count 파일 다운로드
getGEOSuppFiles("GSE255484", makeDirectory = TRUE)

# 3) raw count 읽기
count_df <- read.delim(
  "GSE255484/GSE255484_raw_counts.txt.gz",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(count_df)
mode(expr) <- "numeric"

# 4) expr 열이름에서 내부 샘플 ID 추출
# 예: "4_166__EPG" -> "166"
expr_id <- sub("^[0-9]+_([0-9]+)__EPG$", "\\1", colnames(expr))

# 5) meta$title 에서 내부 샘플 ID 추출
# 예: "... id 166" -> "166"
meta$id_in_title <- sub(".*id[[:space:]]+([0-9]+).*", "\\1", meta$title)

# 6) expr 순서에 맞게 meta 재정렬
meta <- meta[match(expr_id, meta$id_in_title), ]
rownames(meta) <- colnames(expr)

# 7) 샘플 매칭 확인
cat("샘플 이름 일치 여부: ", all(colnames(expr) == rownames(meta)), "\n")
cat("매칭 실패 개수: ", sum(is.na(meta$id_in_title)), "\n")

# 8) 그룹 생성
# responder -> Responder
# control -> Control
# stable + non-responder -> Other
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

# 9) 기본 확인
cat("expr 차원: ", dim(expr)[1], "x", dim(expr)[2], "\n")
cat("meta 차원: ", dim(meta)[1], "x", dim(meta)[2], "\n")
cat("양수 count 개수: ", sum(expr > 0, na.rm = TRUE), "\n")
print(table(meta$group))

# 10) DESeq2 객체 생성
dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = meta,
  design = ~ group
)

# 11) 낮은 발현 유전자 제거
dds <- dds[rowSums(counts(dds)) > 10, ]

# 12) DESeq 실행
dds <- DESeq(dds)

# 13) 결과 추출
res_RC <- results(dds, contrast = c("group", "Responder", "Control"))
res_OC <- results(dds, contrast = c("group", "Other", "Control"))
res_RO <- results(dds, contrast = c("group", "Responder", "Other"))

# 14) 결과 정렬
res_RC_df <- as.data.frame(res_RC[order(res_RC$padj), ])
res_OC_df <- as.data.frame(res_OC[order(res_OC$padj), ])
res_RO_df <- as.data.frame(res_RO[order(res_RO$padj), ])

# gene 이름 컬럼 추가
res_RC_df$gene <- rownames(res_RC_df)
res_OC_df$gene <- rownames(res_OC_df)
res_RO_df$gene <- rownames(res_RO_df)

# 15) 유의한 유전자 추출
sig_RC <- subset(res_RC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_OC <- subset(res_OC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_RO <- subset(res_RO_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# 16) 결과 저장
write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(res_OC_df, "DEG_Other_vs_Control.csv", row.names = FALSE)
write.csv(res_RO_df, "DEG_Responder_vs_Other.csv", row.names = FALSE)

write.csv(sig_RC, "SIG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(sig_OC, "SIG_Other_vs_Control.csv", row.names = FALSE)
write.csv(sig_RO, "SIG_Responder_vs_Other.csv", row.names = FALSE)

# ============================================
# 17) PCA
# ============================================

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA plot") +
  theme_bw()

ggsave("PCA_plot.png", plot = p_pca, width = 8, height = 6, dpi = 300)

# ============================================
# 18) Heatmap
# ============================================

# 분산이 큰 상위 50개 유전자 사용
mat <- assay(vsd)
gene_var <- apply(mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, nrow(mat))]
mat_top <- mat[top_genes, ]

# 유전자별 z-score 변환
mat_top_scaled <- t(scale(t(mat_top)))

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- rownames(meta)

png("Heatmap_top50_genes.png", width = 1200, height = 1000, res = 150)
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

# ============================================
# 19) Volcano plot 함수
# ============================================

make_volcano <- function(res_df, title_text, out_file, top_n = 10) {
  df <- res_df
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  df$significance <- "NS"
  df$significance[df$padj < 0.05 & df$log2FoldChange >= 1] <- "Up"
  df$significance[df$padj < 0.05 & df$log2FoldChange <= -1] <- "Down"
  
  df$neglog10padj <- -log10(df$padj)
  
  # 라벨용 top gene 선택
  df_label <- df[df$padj < 0.05 & abs(df$log2FoldChange) >= 1, ]
  df_label <- df_label[order(df_label$padj), ]
  df_label <- head(df_label, top_n)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj, color = significance)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggtitle(title_text) +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted p-value") +
    theme_bw()
  
  if (nrow(df_label) > 0 && has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(
      data = df_label,
      aes(label = gene),
      size = 3,
      show.legend = FALSE,
      max.overlaps = 20
    )
  }
  
  ggsave(out_file, plot = p, width = 8, height = 6, dpi = 300)
}

# 20) Volcano plot 생성
make_volcano(
  res_RC_df,
  "Volcano: Responder vs Control",
  "Volcano_Responder_vs_Control.png"
)

make_volcano(
  res_OC_df,
  "Volcano: Other vs Control",
  "Volcano_Other_vs_Control.png"
)

make_volcano(
  res_RO_df,
  "Volcano: Responder vs Other",
  "Volcano_Responder_vs_Other.png"
)

# ============================================
# 21) 상위 DEG heatmap (Responder vs Control 기준)
# ============================================

top_sig_RC <- sig_RC[order(sig_RC$padj), ]
top_sig_genes <- head(top_sig_RC$gene, 30)

if (length(top_sig_genes) >= 2) {
  mat_sig <- assay(vsd)[top_sig_genes, , drop = FALSE]
  mat_sig_scaled <- t(scale(t(mat_sig)))
  
  png("Heatmap_top30_SIG_Responder_vs_Control.png", width = 1200, height = 1000, res = 150)
  pheatmap(
    mat_sig_scaled,
    annotation_col = annotation_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 8,
    main = "Top 30 significant genes: Responder vs Control"
  )
  dev.off()
}

# ============================================
# 22) 분석 완료 메시지
# ============================================

cat("\n===== 분석 완료 =====\n")
cat("필터링 후 dds 차원: ", dim(dds)[1], "x", dim(dds)[2], "\n")
cat("Responder vs Control 유의 유전자 수: ", nrow(sig_RC), "\n")
cat("Other vs Control 유의 유전자 수: ", nrow(sig_OC), "\n")
cat("Responder vs Other 유의 유전자 수: ", nrow(sig_RO), "\n")

cat("\n저장된 파일:\n")
cat("- DEG_Responder_vs_Control.csv\n")
cat("- DEG_Other_vs_Control.csv\n")
cat("- DEG_Responder_vs_Other.csv\n")
cat("- SIG_Responder_vs_Control.csv\n")
cat("- SIG_Other_vs_Control.csv\n")
cat("- SIG_Responder_vs_Other.csv\n")
cat("- PCA_plot.png\n")
cat("- Heatmap_top50_genes.png\n")
cat("- Heatmap_top30_SIG_Responder_vs_Control.png\n")
cat("- Volcano_Responder_vs_Control.png\n")
cat("- Volcano_Other_vs_Control.png\n")
cat("- Volcano_Responder_vs_Other.png\n")