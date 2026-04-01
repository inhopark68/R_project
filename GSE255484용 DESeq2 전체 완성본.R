# =========================
# GSE255484 DESeq2 전체 분석 코드
# =========================

# 0) 패키지 로드
library(GEOquery)
library(DESeq2)

# 1) GEO 메타데이터 불러오기
gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

# 확인
head(meta$title, 10)
dim(meta)

# 2) raw counts 파일 다운로드
# 이미 받아졌으면 다시 안 받아도 됨
getGEOSuppFiles("GSE255484", makeDirectory = TRUE)

# 3) raw counts 읽기
count_df <- read.delim(
  "GSE255484/GSE255484_raw_counts.txt.gz",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(count_df)
mode(expr) <- "numeric"

# 확인
dim(expr)
head(expr[, 1:min(5, ncol(expr))])

# 4) expr 열이름에서 내부 샘플 ID 추출
# 예: "4_166__EPG" -> "166"
expr_id <- sub("^[0-9]+_([0-9]+)__EPG$", "\\1", colnames(expr))

# 5) meta$title에서 내부 샘플 ID 추출
# 예: "tumor, anti-CTLA-4, responder, id 166" -> "166"
meta$id_in_title <- sub(".*id[[:space:]]+([0-9]+).*", "\\1", meta$title)

# 6) expr 순서에 맞게 meta 재정렬
meta2 <- meta[match(expr_id, meta$id_in_title), ]
rownames(meta2) <- colnames(expr)

# 매칭 확인
print(all(expr_id == meta2$id_in_title))
print(sum(is.na(meta2$id_in_title)))

# 7) meta 갱신
meta <- meta2

# 8) 그룹 생성
# responder = Responder
# control = Control
# stable + non-responder = Other
meta$group <- ifelse(
  grepl("responder", meta$title, ignore.case = TRUE), "Responder",
  ifelse(
    grepl("control", meta$title, ignore.case = TRUE), "Control", "Other"
  )
)

meta$group <- factor(meta$group, levels = c("Control", "Other", "Responder"))

# 확인
print(all(colnames(expr) == rownames(meta)))
print(table(meta$group))

# 9) DESeq2 객체 생성
dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = meta,
  design = ~ group
)

# 10) 낮은 발현 유전자 제거
dds <- dds[rowSums(counts(dds)) > 10, ]

# 11) DESeq 실행
dds <- DESeq(dds)

# 12) 비교 결과 추출
# Responder vs Control
res_RC <- results(dds, contrast = c("group", "Responder", "Control"))

# Other vs Control
res_OC <- results(dds, contrast = c("group", "Other", "Control"))

# Responder vs Other
res_RO <- results(dds, contrast = c("group", "Responder", "Other"))

# 13) 정렬해서 상위 결과 보기
res_RC_df <- as.data.frame(res_RC[order(res_RC$padj), ])
res_OC_df <- as.data.frame(res_OC[order(res_OC$padj), ])
res_RO_df <- as.data.frame(res_RO[order(res_RO$padj), ])

head(res_RC_df)
head(res_OC_df)
head(res_RO_df)

# 14) 유의한 유전자 추출 예시
sig_RC <- subset(res_RC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_OC <- subset(res_OC_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
sig_RO <- subset(res_RO_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

nrow(sig_RC)
nrow(sig_OC)
nrow(sig_RO)

# 15) 결과 저장
write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = TRUE)
write.csv(res_OC_df, "DEG_Other_vs_Control.csv", row.names = TRUE)
write.csv(res_RO_df, "DEG_Responder_vs_Other.csv", row.names = TRUE)

write.csv(sig_RC, "SIG_Responder_vs_Control.csv", row.names = TRUE)
write.csv(sig_OC, "SIG_Other_vs_Control.csv", row.names = TRUE)
write.csv(sig_RO, "SIG_Responder_vs_Other.csv", row.names = TRUE)

# 16) 간단한 점검 출력
cat("expr dim:", dim(expr), "\n")
cat("meta dim:", dim(meta), "\n")
cat("dds dim:", dim(dds), "\n")
cat("Positive counts:", sum(expr > 0, na.rm = TRUE), "\n")
cat("Groups:\n")
print(table(meta$group))