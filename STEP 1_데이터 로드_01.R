library(GEOquery)

gse <- getGEO("GSE255484", GSEMatrix = TRUE)
eset <- gse[[1]]
meta <- pData(eset)

getGEOSuppFiles("GSE255484", makeDirectory = TRUE)

count_df <- read.delim(
  "GSE255484/GSE255484_raw_counts.txt.gz",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(count_df)
mode(expr) <- "numeric"