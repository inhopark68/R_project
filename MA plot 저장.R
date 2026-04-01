# =========================
# MA plot 저장
# =========================

png("MAplot_Responder_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(
  results(dds, contrast = c("group", "Responder", "Control")),
  main = "MA plot: Responder vs Control",
  ylim = c(-5, 5)
)
dev.off()

png("MAplot_Other_vs_Control.png", width = 1200, height = 1000, res = 150)
plotMA(
  results(dds, contrast = c("group", "Other", "Control")),
  main = "MA plot: Other vs Control",
  ylim = c(-5, 5)
)
dev.off()

png("MAplot_Responder_vs_Other.png", width = 1200, height = 1000, res = 150)
plotMA(
  results(dds, contrast = c("group", "Responder", "Other")),
  main = "MA plot: Responder vs Other",
  ylim = c(-5, 5)
)
dev.off()