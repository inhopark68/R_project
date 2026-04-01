# ============================================
# 관심 유전자 DEG 결과 추출
# ============================================

focus_RC2 <- res_RC_df[res_RC_df$gene_symbol %in% genes_combined, ]
focus_OC2 <- res_OC_df[res_OC_df$gene_symbol %in% genes_combined, ]
focus_RO2 <- res_RO_df[res_RO_df$gene_symbol %in% genes_combined, ]

focus_RC2 <- focus_RC2[order(focus_RC2$padj), ]
focus_OC2 <- focus_OC2[order(focus_OC2$padj), ]
focus_RO2 <- focus_RO2[order(focus_RO2$padj), ]

write.csv(focus_RC2, "FOCUS_Immune_Melanogenesis_Responder_vs_Control.csv", row.names = FALSE)
write.csv(focus_OC2, "FOCUS_Immune_Melanogenesis_Other_vs_Control.csv", row.names = FALSE)
write.csv(focus_RO2, "FOCUS_Immune_Melanogenesis_Responder_vs_Other.csv", row.names = FALSE)

focus_RO2[, c("gene_symbol", "log2FoldChange", "padj")]