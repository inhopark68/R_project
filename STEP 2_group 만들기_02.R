expr_id <- sub("^[0-9]+_([0-9]+)__EPG$", "\\1", colnames(expr))
meta$id_in_title <- sub(".*id[[:space:]]+([0-9]+).*", "\\1", meta$title)

meta <- meta[match(expr_id, meta$id_in_title), ]
rownames(meta) <- colnames(expr)

meta$group <- ifelse(
  grepl("responder", meta$title, ignore.case = TRUE), "Responder",
  ifelse(
    grepl("control", meta$title, ignore.case = TRUE), "Control", "Other"
  )
)

meta$group <- factor(meta$group, levels = c("Control", "Other", "Responder"))