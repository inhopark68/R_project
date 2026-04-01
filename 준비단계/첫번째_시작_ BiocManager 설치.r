# 1) BiocManager 설치
install.packages("BiocManager")

# 2) Bioconductor 버전 확인
BiocManager::version()

# 3) GEOquery와 의존 패키지 설치
BiocManager::install(c("GEOquery", "Biobase"), ask = FALSE, update = TRUE)

# 4) 로드
library(GEOquery)
library(Biobase)