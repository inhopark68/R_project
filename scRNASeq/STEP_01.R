library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

# =========================================================
# 0) input assumptions
# =========================================================
# seu: Seurat object
# metadata columns assumed:
#   - celltype : contains "Malignant"
#   - treatment: e.g. "IgG", "PD1", "CTLA4", "PD1+CTLA4", etc.
#
# If your metadata names differ, change them below.

# =========================================================
# 1) subset malignant cells
# =========================================================
mal <- subset(seu, subset = celltype == "Malignant")

cat("Number of malignant cells:", ncol(mal), "\n")
print(table(mal$treatment))

# =========================================================
# 2) gene sets (mouse symbols)
# =========================================================
melanogenesis_genes <- c("Mitf", "Tyr", "Dct", "Pmel", "Mlana")
dediff_genes       <- c("Axl", "Ngfr", "Egfr", "Sox9", "Zeb1")
ifn_genes          <- c("Stat1", "Ifit1", "Ifit3", "Isg15", "Cxcl10")

# optional extended sets for sensitivity check
melanogenesis_genes_ext <- c("Mitf", "Tyr", "Dct", "Pmel", "Mlana", "Tyrp1", "Sox10", "Trpm1", "Ednrb")
dediff_genes_ext       <- c("Axl", "Ngfr", "Egfr", "Sox9", "Zeb1", "Wnt5a", "Jun", "Fosl1", "Fn1")
ifn_genes_ext          <- c("Stat1", "Ifit1", "Ifit2", "Ifit3", "Isg15", "Ifi6", "Cxcl9", "Cxcl10", "Cxcl11", "B2m")

# keep only genes present in the object
melanogenesis_genes <- intersect(melanogenesis_genes, rownames(mal))
dediff_genes        <- intersect(dediff_genes, rownames(mal))
ifn_genes           <- intersect(ifn_genes, rownames(mal))

cat("Melanogenesis genes used:\n")
print(melanogenesis_genes)

cat("Dediff genes used:\n")
print(dediff_genes)

cat("IFN genes used:\n")
print(ifn_genes)

# safety check
stopifnot(length(melanogenesis_genes) >= 3)
stopifnot(length(dediff_genes) >= 3)
stopifnot(length(ifn_genes) >= 3)

# =========================================================
# 3) module scores
# =========================================================
mal <- AddModuleScore(
  object = mal,
  features = list(melanogenesis_genes),
  name = "melanogenesis_score"
)

mal <- AddModuleScore(
  object = mal,
  features = list(dediff_genes),
  name = "dediff_score"
)

mal <- AddModuleScore(
  object = mal,
  features = list(ifn_genes),
  name = "ifn_score"
)

# AddModuleScore creates columns ending with "1"
score_cols <- c("melanogenesis_score1", "dediff_score1", "ifn_score1")
print(score_cols)

# =========================================================
# 4) z-score normalization across malignant cells
# =========================================================
meta <- mal@meta.data

meta$mel_z    <- as.numeric(scale(meta$melanogenesis_score1))
meta$dediff_z <- as.numeric(scale(meta$dediff_score1))
meta$ifn_z    <- as.numeric(scale(meta$ifn_score1))

# =========================================================
# 5) rule-based state assignment
# =========================================================
meta$cell_state <- "Intermediate"

meta$cell_state[
  meta$mel_z > 0.5 & meta$dediff_z < 0
] <- "MITF_high"

meta$cell_state[
  meta$dediff_z > 0.5 & meta$ifn_z > 0 & meta$mel_z < 0
] <- "AXL_high_inflammatory"

meta$cell_state <- factor(
  meta$cell_state,
  levels = c("MITF_high", "Intermediate", "AXL_high_inflammatory")
)

mal$cell_state <- meta$cell_state
mal$mel_z <- meta$mel_z
mal$dediff_z <- meta$dediff_z
mal$ifn_z <- meta$ifn_z

print(table(mal$cell_state))
print(prop.table(table(mal$cell_state)))

# =========================================================
# 6) state proportion by treatment
# =========================================================
prop_df <- mal@meta.data %>%
  count(treatment, cell_state) %>%
  group_by(treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# fill missing combinations
prop_df <- prop_df %>%
  complete(
    treatment,
    cell_state = factor(c("MITF_high", "Intermediate", "AXL_high_inflammatory"),
                        levels = c("MITF_high", "Intermediate", "AXL_high_inflammatory")),
    fill = list(n = 0, prop = 0)
  )

# order treatments if desired
prop_df$treatment <- factor(
  prop_df$treatment,
  levels = c("IgG", "PD1", "CTLA4", "PD1+Lag3", "PD1+CTLA4", "Lag3")
)

p_bar <- ggplot(prop_df, aes(x = treatment, y = prop, fill = cell_state)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Proportion of malignant cells",
    fill = "Cell state"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar)

# =========================================================
# 7) score summaries by treatment
# =========================================================
score_summary <- mal@meta.data %>%
  group_by(treatment) %>%
  summarise(
    n_cells = n(),
    melanogenesis_mean = mean(melanogenesis_score1, na.rm = TRUE),
    dediff_mean        = mean(dediff_score1, na.rm = TRUE),
    ifn_mean           = mean(ifn_score1, na.rm = TRUE),
    mel_z_mean         = mean(mel_z, na.rm = TRUE),
    dediff_z_mean      = mean(dediff_z, na.rm = TRUE),
    ifn_z_mean         = mean(ifn_z, na.rm = TRUE)
  )

print(score_summary)

# =========================================================
# 8) UMAP panels
# =========================================================
p_umap_state <- DimPlot(
  mal,
  group.by = "cell_state",
  reduction = "umap",
  pt.size = 0.4
) + ggtitle("Malignant cell states")

p_umap_treat <- DimPlot(
  mal,
  group.by = "treatment",
  reduction = "umap",
  pt.size = 0.4
) + ggtitle("Treatment")

print(p_umap_state)
print(p_umap_treat)

# score feature plots
FeaturePlot(
  mal,
  features = c("melanogenesis_score1", "dediff_score1", "ifn_score1"),
  reduction = "umap",
  pt.size = 0.4
)

# =========================================================
# 9) violin plots (supporting only)
# =========================================================
VlnPlot(
  mal,
  features = c("melanogenesis_score1", "dediff_score1", "ifn_score1"),
  group.by = "treatment",
  pt.size = 0
)

# =========================================================
# 10) export state proportion table
# =========================================================
write.csv(prop_df, "GSE243281_malignant_state_proportion.csv", row.names = FALSE)
write.csv(score_summary, "GSE243281_malignant_score_summary.csv", row.names = FALSE)

# =========================================================
# 11) optional: save annotated object
# =========================================================
saveRDS(mal, file = "GSE243281_malignant_with_states.rds")