# ============================================
# GSE255484 DESeq2 전체 분석 최종 스크립트
# ============================================

# 0) 패키지 로드
library(GEOquery)
library(DESeq2)

# 1) GEO 메타데이터 불러오기
gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

# 2) raw count 파일 다운로드
# 이미 있으면 다시 다운로드하지 않아도 됩니다.
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

# 15) 유의한 유전자 추출
sig_RC <- subset(res_RC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_OC <- subset(res_OC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_RO <- subset(res_RO_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# 16) 결과 저장
write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = TRUE)
write.csv(res_OC_df, "DEG_Other_vs_Control.csv", row.names = TRUE)
write.csv(res_RO_df, "DEG_Responder_vs_Other.csv", row.names = TRUE)

write.csv(sig_RC, "SIG_Responder_vs_Control.csv", row.names = TRUE)
write.csv(sig_OC, "SIG_Other_vs_Control.csv", row.names = TRUE)
write.csv(sig_RO, "SIG_Responder_vs_Other.csv", row.names = TRUE)

# 17) 결과 확인
cat("\n===== 분석 완료 =====\n")
cat("필터링 후 dds 차원: ", dim(dds)[1], "x", dim(dds)[2], "\n")
cat("Responder vs Control 유의 유전자 수: ", nrow(sig_RC), "\n")
cat("Other vs Control 유의 유전자 수: ", nrow(sig_OC), "\n")
cat("Responder vs Other 유의 유전자 수: ", nrow(sig_RO), "\n")

cat("\n상위 결과 예시: Responder vs Control\n")
print(head(res_RC_df))

cat("\n저장된 파일:\n")
cat("- DEG_Responder_vs_Control.csv\n")
cat("- DEG_Other_vs_Control.csv\n")
cat("- DEG_Responder_vs_Other.csv\n")
cat("- SIG_Responder_vs_Control.csv\n")
cat("- SIG_Other_vs_Control.csv\n")
cat("- SIG_Responder_vs_Other.csv\n")