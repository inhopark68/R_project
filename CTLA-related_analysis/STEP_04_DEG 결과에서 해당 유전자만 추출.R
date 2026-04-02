# ============================================
# 관심 유전자 DEG 결과만 추출
# ============================================

focus_RC <- res_RC_df[res_RC_df$gene_symbol %in% genes_of_interest, ]
focus_OC <- res_OC_df[res_OC_df$gene_symbol %in% genes_of_interest, ]
focus_RO <- res_RO_df[res_RO_df$gene_symbol %in% genes_of_interest, ]

focus_RC <- focus_RC[order(focus_RC$padj), ]
focus_OC <- focus_OC[order(focus_OC$padj), ]
focus_RO <- focus_RO[order(focus_RO$padj), ]

write.csv(focus_RC, "FOCUS_CTLA4_Melanogenesis_Responder_vs_Control.csv", row.names = FALSE)
write.csv(focus_OC, "FOCUS_CTLA4_Melanogenesis_Other_vs_Control.csv", row.names = FALSE)
write.csv(focus_RO, "FOCUS_CTLA4_Melanogenesis_Responder_vs_Other.csv", row.names = FALSE)

focus_RC[, c("gene_symbol", "log2FoldChange", "pvalue", "padj")]
focus_OC[, c("gene_symbol", "log2FoldChange", "pvalue", "padj")]
focus_RO[, c("gene_symbol", "log2FoldChange", "pvalue", "padj")]