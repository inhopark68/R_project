meta$group <- ifelse(grepl("Responder", meta$title, ignore.case = TRUE), "Responder",
              ifelse(grepl("Non", meta$title, ignore.case = TRUE), "NonResponder",
              ifelse(grepl("IgG|control", meta$title, ignore.case = TRUE), "IgG", "Other")))

table(meta$group)