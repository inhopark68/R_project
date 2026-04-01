# gene symbol 대체 포함

# ============================================
# STEP 4_결과 저장_04.R
# DEG 결과 추출 + gene symbol 매핑 + 저장
# ============================================

library(AnnotationDbi)

# -----------------------------
# 1) 종 선택
# -----------------------------
# 이 데이터는 보통 mouse 분석으로 진행합니다.
# 만약 gene symbol이 대부분 NA로 나오면
# org.Hs.eg.db 로 바꿔서 다시 실행해 보세요.

species_db <- "mouse"

if (species_db == "mouse") {
  library(org.Mm.eg.db)
  anno_db <- org.Mm.eg.db
} else {
  library(org.Hs.eg.db)
  anno_db <- org.Hs.eg.db
}

# -----------------------------
# 2) 결과 추출
# -----------------------------
res_RC <- results(dds, contrast = c("group", "Responder", "Control"))
res_OC <- results(dds, contrast = c("group", "Other", "Control"))
res_RO <- results(dds, contrast = c("group", "Responder", "Other"))

# -----------------------------
# 3) gene symbol 변환 함수
# -----------------------------
add_gene_annotation <- function(res_obj, anno_db) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)

  # Ensembl version 제거
  # 예: ENSMUSG00000000001.1 -> ENSMUSG00000000001
  res_df$gene_id_clean <- sub("\\..*$", "", res_df$gene_id)

  # SYMBOL
  res_df$gene_symbol <- mapIds(
    anno_db,
    keys = res_df$gene_id_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  # ENTREZID
  res_df$entrez_id <- mapIds(
    anno_db,
    keys = res_df$gene_id_clean,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  # GENENAME
  res_df$gene_name <- mapIds(
    anno_db,
    keys = res_df$gene_id_clean,
    column = "GENENAME",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  # gene_symbol 없는 경우 원래 gene_id 유지용 보조 컬럼
  res_df$gene_label <- ifelse(
    is.na(res_df$gene_symbol) | res_df$gene_symbol == "",
    res_df$gene_id,
    res_df$gene_symbol
  )

  # padj 기준 정렬
  res_df <- res_df[order(res_df$padj), ]

  return(res_df)
}

# -----------------------------
# 4) annotation 추가
# -----------------------------
res_RC_df <- add_gene_annotation(res_RC, anno_db)
res_OC_df <- add_gene_annotation(res_OC, anno_db)
res_RO_df <- add_gene_annotation(res_RO, anno_db)

# -----------------------------
# 5) 유의 유전자 추출
# -----------------------------
sig_RC <- subset(res_RC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_OC <- subset(res_OC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_RO <- subset(res_RO_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# -----------------------------
# 6) 저장
# -----------------------------
write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(res_OC_df, "DEG_Other_vs_Control.csv", row.names = FALSE)
write.csv(res_RO_df, "DEG_Responder_vs_Other.csv", row.names = FALSE)

write.csv(sig_RC, "SIG_Responder_vs_Control.csv", row.names = FALSE)
write.csv(sig_OC, "SIG_Other_vs_Control.csv", row.names = FALSE)
write.csv(sig_RO, "SIG_Responder_vs_Other.csv", row.names = FALSE)

# -----------------------------
# 7) 확인
# -----------------------------
cat("Responder vs Control significant:", nrow(sig_RC), "\n")
cat("Other vs Control significant:", nrow(sig_OC), "\n")
cat("Responder vs Other significant:", nrow(sig_RO), "\n")

cat("\nExample annotation:\n")
print(head(res_RC_df[, c("gene_id", "gene_symbol", "gene_name", "log2FoldChange", "padj")]))