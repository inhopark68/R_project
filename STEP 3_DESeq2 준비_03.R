library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = meta,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)