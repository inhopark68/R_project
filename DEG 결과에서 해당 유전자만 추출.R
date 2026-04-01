# Responder vs Control 기준
subset_RC <- res_RC_df[res_RC_df$gene_symbol %in% genes_of_interest, ]

# Other vs Control
subset_OC <- res_OC_df[res_OC_df$gene_symbol %in% genes_of_interest, ]

# Responder vs Other
subset_RO <- res_RO_df[res_RO_df$gene_symbol %in% genes_of_interest, ]

subset_RC