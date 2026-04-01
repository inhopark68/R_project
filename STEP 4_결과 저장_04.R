res_RC <- results(dds, contrast = c("group", "Responder", "Control"))
res_OC <- results(dds, contrast = c("group", "Other", "Control"))
res_RO <- results(dds, contrast = c("group", "Responder", "Other"))

res_RC_df <- as.data.frame(res_RC[order(res_RC$padj), ])
res_OC_df <- as.data.frame(res_OC[order(res_OC$padj), ])
res_RO_df <- as.data.frame(res_RO[order(res_RO$padj), ])

write.csv(res_RC_df, "DEG_Responder_vs_Control.csv", row.names = TRUE)
write.csv(res_OC_df, "DEG_Other_vs_Control.csv", row.names = TRUE)
write.csv(res_RO_df, "DEG_Responder_vs_Other.csv", row.names = TRUE)