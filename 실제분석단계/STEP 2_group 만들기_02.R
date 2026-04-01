meta$group <- ifelse(grepl("Responder", meta$title, ignore.case = TRUE), "Responder",
              ifelse(grepl("Non", meta$title, ignore.case = TRUE), "NonResponder",
              ifelse(grepl("IgG|control", meta$title, ignore.case = TRUE), "Control", "Other")))

table(meta$group)

#   Control     Other Responder 
#         5         6        14 

# 👉 여기서 반드시 확인:

# Responder 몇 개
# NonResponder 몇 개
# Control 몇 개