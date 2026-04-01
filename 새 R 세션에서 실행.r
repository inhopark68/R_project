# 1) 새 R 세션에서 실행
remove.packages("SparseArray")

# 2) Bioconductor 버전 확인
BiocManager::version()

# 3) 문제 패키지 다시 설치
BiocManager::install("SparseArray", force = TRUE)

# 4) 그 다음 DESeq2 재확인
BiocManager::install("DESeq2", force = TRUE)

# 5) 로드 테스트
library(DESeq2)