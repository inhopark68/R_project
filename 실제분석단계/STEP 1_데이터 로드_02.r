eset <- gse[[1]]

expr <- exprs(eset)
meta <- pData(eset)

dim(expr)

# [1]  0 25

dim(meta)

# [1] 25 46

head(meta[,1:10])

                                        #    title geo_accession
# GSM8073211 tumor, anti-CTLA-4, responder, id 140    GSM8073211
# GSM8073212 tumor, anti-CTLA-4, responder, id 142    GSM8073212
# GSM8073213 tumor, anti-CTLA-4, responder, id 149    GSM8073213
# GSM8073214 tumor, anti-CTLA-4, responder, id 160    GSM8073214
# GSM8073215 tumor, anti-CTLA-4, responder, id 161    GSM8073215
# GSM8073216 tumor, anti-CTLA-4, responder, id 163    GSM8073216
#                           status submission_date last_update_date type
# GSM8073211 Public on Feb 29 2024     Feb 09 2024      Feb 29 2024  SRA
# GSM8073212 Public on Feb 29 2024     Feb 09 2024      Feb 29 2024  SRA
# GSM8073213 Public on Feb 29 2024     Feb 09 2024      Feb 29 2024  SRA
# GSM8073214 Public on Feb 29 2024     Feb 09 2024      Feb 29 2024  SRA
# GSM8073215 Public on Feb 29 2024     Feb 09 2024      Feb 29 2024  SRA
# GSM8073216 Public on Feb 29 2024     Feb 09 2024      Feb 29 2024  SRA
#            channel_count source_name_ch1 organism_ch1 characteristics_ch1
# GSM8073211             1           B2905 Mus musculus    cell line: B2905
# GSM8073212             1           B2905 Mus musculus    cell line: B2905
# GSM8073213             1           B2905 Mus musculus    cell line: B2905
# GSM8073214             1           B2905 Mus musculus    cell line: B2905
# GSM8073215             1           B2905 Mus musculus    cell line: B2905
# GSM8073216             1           B2905 Mus musculus    cell line: B2905