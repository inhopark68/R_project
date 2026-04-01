library(GEOquery)
library(Biobase)

gse <- getGEO("GSE255484", GSEMatrix = TRUE)

length(gse)

eset <- gse[[1]]
expr <- exprs(eset)
meta <- pData(eset)

dim(expr)
dim(meta)
head(meta[, 1:10])