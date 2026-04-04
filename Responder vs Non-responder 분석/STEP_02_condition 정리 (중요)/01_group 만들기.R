colnames(meta)
unique(meta$title)



meta$group <- ifelse(grepl("Responder", meta$title, ignore.case = TRUE), "Responder",
              ifelse(grepl("Non", meta$title, ignore.case = TRUE), "NonResponder",
              ifelse(grepl("IgG|control", meta$title, ignore.case = TRUE), "Control", "Other")))

table(meta$group)
levels(as.factor(meta$group))
resultsNames(dds)

# 1) 가장 먼저 group 다시 확인
unique(meta$group)
table(meta$group, useNA = "ifany")


# 2) 안전하게 group 다시 만들기
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

# 그다음 Other는 빼고 진행하는 것이 좋습니다.

keep <- meta$group %in% c("Control", "Responder", "NonResponder")
meta2 <- meta[keep, ]
expr2 <- expr[, rownames(meta2)]

# 3) factor level 명시적으로 지정
meta2$group <- factor(meta2$group, levels = c("Control", "NonResponder", "Responder"))
table(meta2$group)
levels(meta2$group)


# 4) DESeq2 다시 만들기
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(expr2),
  colData = meta2,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

resultsNames(dds)


# 5) 이제 contrast 실행
res_R_vs_NR <- results(dds, contrast = c("group", "Responder", "NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)

head(res_R_vs_NR)


# 1. 현재 group 상태 확인

table(meta$group, useNA = "ifany")
unique(meta$group)

# 2. dds 안의 실제 비교 이름 확인
resultsNames(dds)
colData(dds)$group
levels(colData(dds)$group)
table(colData(dds)$group, useNA = "ifany")

# 3. 가장 안전하게 group 다시 만들기
# 기존 코드는 Responder를 먼저 잡아서 NonResponder가 잘못 분류됐을 가능성이 큽니다.
meta$group <- ifelse(
  grepl("NonResponder|Non-responder|non responder|non-responder|NR", meta$title, ignore.case = TRUE),
  "NonResponder",
  ifelse(
    grepl("Responder", meta$title, ignore.case = TRUE),
    "Responder",
    ifelse(
      grepl("IgG|control|Ctrl", meta$title, ignore.case = TRUE),
      "Control",
      "Other"
    )
  )
)

table(meta$group, useNA = "ifany")


# 4. 필요한 샘플만 남기기
keep <- meta$group %in% c("Control", "NonResponder", "Responder")
meta2 <- meta[keep, ]
expr2 <- expr[, rownames(meta2)]

# 5. factor level을 명시적으로 지정
meta2$group <- factor(meta2$group, levels = c("Control", "NonResponder", "Responder"))

table(meta2$group)
levels(meta2$group)

