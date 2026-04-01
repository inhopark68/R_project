focus_RC <- res_RC_df[res_RC_df$gene_symbol %in% genes_of_interest, ]
focus_OC <- res_OC_df[res_OC_df$gene_symbol %in% genes_of_interest, ]
focus_RO <- res_RO_df[res_RO_df$gene_symbol %in% genes_of_interest, ]

focus_RC <- focus_RC[order(focus_RC$padj), ]
focus_OC <- focus_OC[order(focus_OC$padj), ]
focus_RO <- focus_RO[order(focus_RO$padj), ]

focus_RC[, c("gene_symbol", "log2FoldChange", "padj")]
focus_OC[, c("gene_symbol", "log2FoldChange", "padj")]
focus_RO[, c("gene_symbol", "log2FoldChange", "padj")]