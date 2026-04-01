# ============================================
# Immune signature genes
# ============================================

immune_genes <- c(
  "Cd8a", "Cd8b1", "Gzmb", "Prf1", "Ifng",
  "Pdcd1", "Lag3", "Tigit", "Havcr2",
  "Foxp3", "Ctla4", "Tbx21", "Stat1",
  "Cxcl9", "Cxcl10", "Cxcl11"
)

melanogenesis_genes <- c(
  "Mitf", "Tyr", "Tyrp1", "Dct",
  "Pmel", "Mlana", "Sox10", "Mc1r", "Kit", "Gpr143"
)

genes_combined <- unique(c(immune_genes, melanogenesis_genes))
genes_combined