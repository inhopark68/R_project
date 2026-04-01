library(GEOquery)
library(Biobase)

gse <- getGEO("GSE255484", GSEMatrix = TRUE)

# Found 1 file(s)
# GSE255484_series_matrix.txt.gz
# Using locally cached version: C:\Users\inhopark\AppData\Local\Temp\RtmpgVnTg0/GSE255484_series_matrix.txt.gz
# Using locally cached version of GPL13112 found here:


length(gse)

# [1] 1