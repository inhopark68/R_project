table(meta$group, useNA = "ifany")
unique(meta$group)

meta$group <- ifelse(
  grepl("NonResponder|Non-responder|non responder|non-responder|NR", meta$title, ignore.case = TRUE),
  "NonResponder",
  ifelse(
    grepl("Responder|responder", meta$title, ignore.case = TRUE),
    "Responder",
    ifelse(
      grepl("IgG|control|Ctrl", meta$title, ignore.case = TRUE),
      "Control",
      "Other"
    )
  )
)

table(meta$group, useNA = "ifany")

keep <- meta$group %in% c("Control", "NonResponder", "Responder")
meta2 <- meta[keep, ]
expr2 <- expr[, rownames(meta2)]

meta2$group <- factor(meta2$group, levels = c("Control", "NonResponder", "Responder"))

dds <- DESeqDataSetFromMatrix(
  countData = round(expr2),
  colData = meta2,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

resultsNames(dds)

res_R_vs_NR <- results(dds, contrast = c("group", "Responder", "NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)

head(res_R_vs_NR)