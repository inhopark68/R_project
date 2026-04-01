library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = meta,
  design = ~ group
)

# converting counts to integer mode
# Error in DESeqDataSet(se, design = design, ignoreRank) : 
#   all samples have 0 counts for all genes. check the counting script.

dds <- dds[rowSums(counts(dds)) > 10,]

# Error: object 'dds' not found

dds <- DESeq(dds)