# ============================================
# STEP 3_DESeq2 준비_03.R
# DESeq2 실행
# ============================================

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = meta,
  design = ~ group
)

# 낮은 발현 유전자 제거
dds <- dds[rowSums(counts(dds)) > 10, ]

# DESeq 실행
dds <- DESeq(dds)

cat("dds dim after filtering:", dim(dds), "\n")