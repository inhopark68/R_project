res_R_vs_NR <- results(dds, contrast = c("group","Responder","NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)

head(res_R_vs_NR)

2) Responder vs NonResponder 에러 원인

"group_Responder_vs_Control"



# 바로 해결하는 코드
# A. rowname에서 gene symbol만 추출: 오류 나옴

expr_log2 <- expr_log
rownames(expr_log2) <- sub(".*_", "", rownames(expr_log2))

head(rownames(expr_log2), 20)

# 매칭 확인:
intersect(mel_genes, rownames(expr_log2))
intersect(dediff_genes, rownames(expr_log2))
intersect(ifn_genes, rownames(expr_log2))

# B. module score 다시 계산
mel_score <- colMeans(expr_log2[intersect(mel_genes, rownames(expr_log2)), , drop = FALSE])
dediff_score <- colMeans(expr_log2[intersect(dediff_genes, rownames(expr_log2)), , drop = FALSE])
ifn_score <- colMeans(expr_log2[intersect(ifn_genes, rownames(expr_log2)), , drop = FALSE])

module_df <- data.frame(
  sample = colnames(expr_log2),
  melanogenesis = mel_score,
  dediff = dediff_score,
  ifn = ifn_score,
  group = meta$group
)

module_df

# C. group을 다시 정확히 만들기

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


# D. NonResponder가 실제로 title에 어떻게 적혀 있는지 확인

unique(meta$title)



# 특히 NonResponder 샘플 제목을 눈으로 확인해 보세요.

grep("Non|NR|responder|IgG|control", unique(meta$title), value = TRUE, ignore.case = TRUE)


# E. 분석용 데이터 다시 만들기

keep <- meta$group %in% c("Control", "NonResponder", "Responder")
meta2 <- meta[keep, ]
expr2 <- expr[, rownames(meta2)]

meta2$group <- factor(meta2$group, levels = c("Control", "NonResponder", "Responder"))

table(meta2$group)
levels(meta2$group)

# F. DESeq2 다시 실행  오류 발생
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(expr2),
  colData = meta2,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

resultsNames(dds)

# G. 이제 비교 실행
res_R_vs_NR <- results(dds, contrast = c("group", "Responder", "NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)

head(res_R_vs_NR)





# module plot도 다시 그리기

library(ggplot2)

ggplot(module_df, aes(x = group, y = melanogenesis)) +
  geom_boxplot() +
  ggtitle("Melanogenesis score")

ggplot(module_df, aes(x = group, y = dediff)) +
  geom_boxplot() +
  ggtitle("Dedifferentiation score")

ggplot(module_df, aes(x = group, y = ifn)) +
  geom_boxplot() +
  ggtitle("IFN response score")


# 가장 먼저 실행할 4줄

# 이 4줄부터 바로 해보세요.
expr_log2 <- expr_log
rownames(expr_log2) <- sub(".*_", "", rownames(expr_log2))

intersect(mel_genes, rownames(expr_log2))
table(meta$group, useNA = "ifany")
unique(meta$title)

# 지금 다시 해야 하는 것

# 당신이 boxplot을 그릴 때 사용한 module_df는 예전 expr_log 기준으로 만든 것이라서,
# 아마 아직도 NaN 버전일 가능성이 큽니다.

# 즉, expr_log2로 다시 module score를 계산한 뒤 plot을 다시 그려야 합니다.

# 아래를 그대로 실행하세요.

mel_score <- colMeans(expr_log2[intersect(mel_genes, rownames(expr_log2)), , drop = FALSE])
dediff_score <- colMeans(expr_log2[intersect(dediff_genes, rownames(expr_log2)), , drop = FALSE])
ifn_score <- colMeans(expr_log2[intersect(ifn_genes, rownames(expr_log2)), , drop = FALSE])

module_df <- data.frame(
  sample = colnames(expr_log2),
  melanogenesis = mel_score,
  dediff = dediff_score,
  ifn = ifn_score,
  group = meta$group
)

module_df


library(ggplot2)

ggplot(module_df, aes(x = group, y = melanogenesis)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  ggtitle("Melanogenesis score")

ggplot(module_df, aes(x = group, y = dediff)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  ggtitle("Dedifferentiation score")

ggplot(module_df, aes(x = group, y = ifn)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  ggtitle("IFN response score")


#   DE 결과도 gene symbol 붙여서 보기

res_R_vs_NR$symbol <- sub(".*_", "", res_R_vs_NR$gene)
head(res_R_vs_NR[, c("symbol", "log2FoldChange", "padj")])



# 관심 유전자만 보기:

genes_of_interest <- c(
  mel_genes,
  dediff_genes,
  ifn_genes
)

goi_res <- res_R_vs_NR[res_R_vs_NR$symbol %in% genes_of_interest, ]
goi_res <- goi_res[order(goi_res$padj), ]

goi_res[, c("symbol", "log2FoldChange", "padj")]


# volcano plot도 symbol 기준으로 라벨링

res_R_vs_NR$log10padj <- -log10(res_R_vs_NR$padj)

label_genes <- c("Mitf","Tyr","Tyrp1","Dct","Pmel","Mlana","Axl","Ngfr","Stat1","Irf1","Ifit1","Isg15")

ggplot(res_R_vs_NR, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_text(
    data = subset(res_R_vs_NR, symbol %in% label_genes),
    aes(label = symbol),
    size = 3,
    vjust = -0.5
  ) +
  theme_classic()



# 지금부터 해석하는 법

"Responder vs NonResponder"

# 2) 실제 비교 결과를 보고 싶으면
head(res_R_vs_NR)

# 3) 관심 유전자만 보고 싶으면

res_R_vs_NR$symbol <- sub(".*_", "", res_R_vs_NR$gene)

genes_of_interest <- c(
  mel_genes,
  dediff_genes,
  ifn_genes
)

goi_res <- res_R_vs_NR[res_R_vs_NR$symbol %in% genes_of_interest, ]
goi_res <- goi_res[order(goi_res$padj), ]

goi_res[, c("symbol", "log2FoldChange", "padj")]

# 지금부터는 이렇게 생각하시면 됩니다

# 당신의 DE 비교는 이미 이 줄로 끝났습니다.

res_R_vs_NR <- results(dds, contrast = c("group", "Responder", "NonResponder"))

# 바로 다음으로 실행할 것
res_R_vs_NR$symbol <- sub(".*_", "", res_R_vs_NR$gene)

genes_of_interest <- c(
  mel_genes,
  dediff_genes,
  ifn_genes
)

goi_res <- res_R_vs_NR[res_R_vs_NR$symbol %in% genes_of_interest, ]
goi_res <- goi_res[order(goi_res$padj), ]

goi_res[, c("symbol", "log2FoldChange", "padj")]



# 바로 해결 코드

# 아래를 그대로 순서대로 실행하세요.

res_R_vs_NR <- results(dds, contrast = c("group", "Responder", "NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)
res_R_vs_NR$symbol <- sub(".*_", "", res_R_vs_NR$gene)

head(res_R_vs_NR[, c("gene", "symbol", "log2FoldChange", "padj")])


# 관심 유전자만 보기

genes_of_interest <- c(
  mel_genes,
  dediff_genes,
  ifn_genes
)

goi_res <- res_R_vs_NR[res_R_vs_NR$symbol %in% genes_of_interest, ]
goi_res <- goi_res[order(goi_res$padj), ]

goi_res[, c("symbol", "log2FoldChange", "padj")]

# 안전하게 한 번에 실행하는 버전

res_R_vs_NR <- results(dds, contrast = c("group", "Responder", "NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)
res_R_vs_NR$symbol <- sub(".*_", "", res_R_vs_NR$gene)

genes_of_interest <- c(
  mel_genes,
  dediff_genes,
  ifn_genes
)

goi_res <- res_R_vs_NR[res_R_vs_NR$symbol %in% genes_of_interest, , drop = FALSE]
goi_res <- goi_res[order(goi_res$padj), , drop = FALSE]

print(goi_res[, c("symbol", "log2FoldChange", "padj"), drop = FALSE])


# 추가로 권장

# 지금 group에 Stable을 따로 두는 게 더 해석에 좋습니다.
# 하지만 우선은 지금처럼 Responder vs NonResponder만 먼저 보는 게 맞습니다.

# 지금 다음으로 필요한 건 이 출력입니다:

goi_res[, c("symbol", "log2FoldChange", "padj")]

# 그 결과를 붙여주시면 제가 바로
# melanogenesis 감소/증가, dediff 증가/감소, IFN 변화를 해석해드리겠습니다.