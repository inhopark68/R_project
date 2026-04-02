res_R_vs_NR <- results(dds, contrast = c("group","Responder","NonResponder"))
res_R_vs_NR <- as.data.frame(res_R_vs_NR)
res_R_vs_NR$gene <- rownames(res_R_vs_NR)

head(res_R_vs_NR)