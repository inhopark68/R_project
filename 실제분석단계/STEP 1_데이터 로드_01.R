# ============================================
# STEP 1_데이터 로드_01.R
# GEO metadata + raw counts 불러오기
# ============================================

library(GEOquery)

# 1) GEO 메타데이터 불러오기
gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

# 2) raw counts 다운로드
getGEOSuppFiles("GSE255484", makeDirectory = TRUE)

# 3) raw counts 읽기
count_df <- read.delim(
  "GSE255484/GSE255484_raw_counts.txt.gz",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(count_df)
mode(expr) <- "numeric"

# 4) 확인
cat("expr dim:", dim(expr), "\n")
cat("meta dim:", dim(meta), "\n")
head(meta$title, 10)
head(colnames(expr), 10)