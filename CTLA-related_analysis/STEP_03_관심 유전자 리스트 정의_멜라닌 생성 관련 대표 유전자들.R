# ============================================
# 관심 유전자 목록
# ============================================

ctla4_gene <- "Ctla4"

melanogenesis_genes <- c(
  "Mitf", "Tyr", "Tyrp1", "Dct",
  "Pmel", "Mlana", "Sox10",
  "Mc1r", "Kit", "Gpr143"
)

genes_of_interest <- unique(c(ctla4_gene, melanogenesis_genes))
genes_of_interest