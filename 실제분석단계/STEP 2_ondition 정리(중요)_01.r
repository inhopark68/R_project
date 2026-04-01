#이 데이터는 title이나 characteristics에 condition이 들어있습니다.
#먼저 확인:

colnames(meta)

#  [1] "title"                   "geo_accession"          
#  [3] "status"                  "submission_date"        
#  [5] "last_update_date"        "type"                   
#  [7] "channel_count"           "source_name_ch1"        
#  [9] "organism_ch1"            "characteristics_ch1"    
# [11] "characteristics_ch1.1"   "treatment_protocol_ch1" 
# [13] "growth_protocol_ch1"     "molecule_ch1"           
# [15] "extract_protocol_ch1"    "extract_protocol_ch1.1" 
# [17] "taxid_ch1"               "description"            
# [19] "data_processing"         "data_processing.1"      
# [21] "data_processing.2"       "data_processing.3"      
# [23] "data_processing.4"       "data_processing.5"      
# [25] "data_processing.6"       "platform_id"            
# [27] "contact_name"            "contact_email"          
# [29] "contact_laboratory"      "contact_department"     
# [31] "contact_institute"       "contact_address"        
# [33] "contact_city"            "contact_state"          
# [35] "contact_zip/postal_code" "contact_country"        
# [37] "data_row_count"          "instrument_model"       
# [39] "library_selection"       "library_source"         
# [41] "library_strategy"        "relation"               
# [43] "relation.1"              "supplementary_file_1"   
# [45] "cell line:ch1"           "treatment:ch1"          
# [47] "group"          

# > 5         6        14
# Error: unexpected numeric constant in "5         6"

# 로그를 보면:

# meta$group → 정상 (Control / Other / Responder 있음) ✅
# meta 구조 → 정상 (group 컬럼 존재) ✅
# 그런데 dds 생성 → 실패 ❌

# 👉 즉, 문제는 100% expr 쪽

# ❌ expr = raw count가 아님

# DESeq2는 반드시 이런 데이터를 요구합니다:

# 정수 count (RNA-seq read count)
# 예: 0, 5, 123, 9876

# 그런데 현재 expr는 아마 이런 상태일 가능성 매우 높습니다:

# log2 값 (예: 6.2, 7.1)
# TPM / FPKM
# normalized 값
# microarray expression 값
# 소수점 데이터

# 👉 그래서 round(expr) 하면 전부 0 → 에러 발생

# 🔥 지금 바로 확인 (중요)

# 이거 그대로 실행해서 결과 보여주세요:

summary(as.matrix(expr))

👉 결과가 0이면 완전히 잘못된 데이터

#  GSM8073211     GSM8073212     GSM8073213     GSM8073214     GSM8073215    
#  Mode:logical   Mode:logical   Mode:logical   Mode:logical   Mode:logical  
#  GSM8073216     GSM8073217     GSM8073218     GSM8073219     GSM8073220    
#  Mode:logical   Mode:logical   Mode:logical   Mode:logical   Mode:logical  
#  GSM8073221     GSM8073222     GSM8073223     GSM8073224     GSM8073225    
#  Mode:logical   Mode:logical   Mode:logical   Mode:logical   Mode:logical  
#  GSM8073226     GSM8073227     GSM8073228     GSM8073229     GSM8073230    
#  Mode:logical   Mode:logical   Mode:logical   Mode:logical   Mode:logical  
#  GSM8073231     GSM8073232     GSM8073233     GSM8073234     GSM8073235    
#  Mode:logical   Mode:logical   Mode:logical   Mode:logical   Mode:logical 

# 원인 찾았습니다.
# 지금 expr는 숫자 행렬이 아니라 logical 행렬입니다.
# 이 뜻은 보통 둘 중 하나입니다.

# expr가 TRUE/FALSE로 만들어졌거나
# expr가 전부 NA가 된 상태입니다

# 그래서 DESeq2가 count 데이터로 쓸 수 없고, all samples have 0 counts for all genes 같은 
# 에러로 이어진 겁니다. 지금까지의 콘솔 로그도 그 흐름과 일치합니다.

# 지금 바로 확인

# 아래를 그대로 실행해 보세요.
class(expr)
[1] "matrix" "array" 
str(expr)
 logi[0 , 1:25] 
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:25] "GSM8073211" "GSM8073212" "GSM8073213" "GSM8073214" ...

dim(expr)
[1]  0 25

head(expr[, 1:min(5, ncol(expr))])
    #  GSM8073211 GSM8073212 GSM8073213 GSM8073214 GSM8073215

table(unlist(expr), useNA = "ifany")
# < table of extent 0 >

# 특히 table(...) 결과가 이런 식이면 상태가 확실해집니다.

# FALSE만 많음 → 논리값 데이터
# TRUE/FALSE 섞임 → 논리 비교 결과를 expr에 넣어버린 것
# NA만 있음 → 변환 과정에서 전부 깨진 것


원인 확정입니다.

지금 expr는 값이 없는 빈 행렬입니다.

class(expr)  # matrix
str(expr)    # logi[0, 1:25]
dim(expr)    # 0 25

즉 뜻이 이겁니다.

열은 25개 샘플이 있음
행은 0개 유전자
타입은 logical
그래서 DESeq2에 넣을 유전자 count가 아예 없음

그래서 DESeqDataSetFromMatrix()가 실패한 겁니다. 
이전 로그의 all samples have 0 counts for all genes와도 일치합니다. 
meta$group은 정상적으로 25개 샘플로 잡혀 있었지만, expr 쪽에 유전자가 하나도 없었습니다.

핵심은 DESeq2 문제가 아니라, expr를 만드는 단계에서 모든 행이 날아갔다는 점입니다.

보통 원인은 셋 중 하나입니다.

1. 필터링을 너무 강하게 해서 전부 제거함
예:
expr <- expr[rowSums(expr > 10) > 5, ]

이런 조건이 너무 강하면 0행이 됩니다.

2. 샘플/행 선택 조건이 전부 FALSE였음
예:
expr <- expr[genes_to_keep, ]

3. 파일을 잘못 읽어서 유전자 행이 처음부터 안 들어옴

지금 바로 확인할 가장 중요한 건 expr를 만드는 직전 객체입니다.
아래를 순서대로 실행해 보세요.

exists("expr")
# [1] TRUE
dim(expr)
# NULL
# expr를 만들기 전 원본 객체 후보 확인
ls()
# [1] "eset" "expr" "gse"  "meta"


# STEP 1 다시 실행
source("d:/R_project/실제분석단계/STEP 1_데이터 로드_01.R", encoding = "UTF-8")

# 실행 직후 확인
class(expr)
dim(expr)
# [1]  0 25
str(expr)
#  logi[0 , 1:25] 
#  - attr(*, "dimnames")=List of 2
#   ..$ : NULL
#   ..$ : chr [1:25] "GSM8073211" "GSM8073212" "GSM8073213" "GSM8073214" ...
head(expr[, 1:min(5, ncol(expr))])
    #  GSM8073211 GSM8073212 GSM8073213 GSM8073214 GSM8073215

여기서 STEP 1을 실행한 직후에도 dim(expr)가 0 25면, 문제는 STEP 1 안에 있습니다.
반대로 STEP 1 직후에는 정상인데 STEP 2 이후 0행이 되면, 문제는 STEP 2 필터링 코드입니다.

가장 빠른 진단법은 각 단계 뒤에 이걸 넣는 겁니다.


이렇게 하면 어느 단계에서 0행이 되는지 바로 잡힙니다.
그리고 DESeq2로 가기 전, 정상 상태의 expr는 최소한 이런 조건이어야 합니다.

is.matrix(expr)
# [1] TRUE
nrow(expr) > 0
# [1] FALSE
ncol(expr) == 25
# [1] TRUE
is.numeric(expr) || is.integer(expr)
# [1] FALSE


원인 찾았습니다.

지금 데이터는 GSE255484이고, GEO 설명에 따르면 “Expression profiling by high throughput sequencing”, 
즉 RNA-seq 연구입니다. 
그리고 이 시리즈에는 별도로 **GSE255484_raw_counts.txt.gz**와 **GSE255484_rpkm_expr.csv.gz**가 제공됩니다.

그래서 지금처럼
eset <- gse[[1]]
expr <- exprs(eset)
를 했을 때 dim(expr) == c(0, 25)가 나온 건, series matrix 안에 실제 expression/count 테이블이 비어 있고 
샘플 정보만 들어왔기 때문입니다. 
당신 콘솔에서도 expr가 logi[0, 1:25]인 빈 행렬로 확인됐습니다.
즉, 문제는 샘플/행 선택이 아니라, 애초에 행(유전자)이 없는 객체를 잡은 것입니다.

1) raw counts 다운로드

library(GEOquery)

getGEOSuppFiles("GSE255484", makeDirectory = TRUE)

그러면 보통 이런 파일이 생깁니다.

GSE255484/GSE255484_raw_counts.txt.gz

#                                                       size isdir mode
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz 1525696 FALSE  666
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  3663023 FALSE  666
#                                                                  mtime
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz 2026-04-01 15:20:05
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  2026-04-01 15:20:14
#                                                                  ctime
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz 2026-04-01 15:19:56
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  2026-04-01 15:20:06
#                                                                  atime exe
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz 2026-04-01 15:20:05  no
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  2026-04-01 15:20:14  no
#                                                       uname udomain
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz INHOPARK  YUHSSC
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  INHOPARK  YUHSSC
#                                                                          fname
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz GSE255484_raw_counts.txt.gz
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz   GSE255484_rpkm_expr.csv.gz
#                                                                   destdir
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz D:/R_project/GSE255484
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  D:/R_project/GSE255484
#                                                                                              filepath
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz   D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz
#                                                          GEO
# D:/R_project/GSE255484/GSE255484_raw_counts.txt.gz GSE255484
# D:/R_project/GSE255484/GSE255484_rpkm_expr.csv.gz  GSE255484

2) count matrix 읽기

count_df <- read.delim(
  "GSE255484/GSE255484_raw_counts.txt.gz",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(count_df)
mode(expr) <- "numeric"

dim(expr)
# [1] 55450    25

head(expr[, 1:5])
#                             4_166__EPG 20_173__EPG 14_149__EPG 18_165__EPG
# ENSMUSG00000000001.4_Gnai3        5383        2217        2403         822
# ENSMUSG00000000003.15_Pbsn           0           0           0           0
# ENSMUSG00000000028.15_Cdc45        637          97         140          78
# ENSMUSG00000000031.16_H19          700         952         329         906
# ENSMUSG00000000037.16_Scml2        212          26          65          38
# ENSMUSG00000000049.11_Apoh           6           6           2           4
#                             5_168__EPG
# ENSMUSG00000000001.4_Gnai3        2920
# ENSMUSG00000000003.15_Pbsn           0
# ENSMUSG00000000028.15_Cdc45        347
# ENSMUSG00000000031.16_H19          488
# ENSMUSG00000000037.16_Scml2         89
# ENSMUSG00000000049.11_Apoh           1


3) meta와 샘플 이름 맞추기

현재 meta는 25개 샘플로 정상입니다. group 분포도 Control 5, Other 6, Responder 14로 잡혀 있었습니다.

이제 샘플 순서만 맞추세요.

meta <- pData(eset)

colnames(expr)
#  [1] "4_166__EPG"  "20_173__EPG" "14_149__EPG" "18_165__EPG" "5_168__EPG" 
#  [6] "16_161__EPG" "2_146__EPG"  "9_150__EPG"  "1_145__EPG"  "21_136__EPG"
# [11] "13_142__EPG" "12_140__EPG" "6_138__EPG"  "22_137__EPG" "8_143__EPG" 
# [16] "15_160__EPG" "19_169__EPG" "25_158__EPG" "7_139__EPG"  "17_163__EPG"
# [21] "24_167__EPG" "23_155__EPG" "11_172__EPG" "10_162__EPG" "3_152__EPG

rownames(meta)
#  [1] "GSM8073211" "GSM8073212" "GSM8073213" "GSM8073214" "GSM8073215"
#  [6] "GSM8073216" "GSM8073217" "GSM8073218" "GSM8073219" "GSM8073220"
# [11] "GSM8073221" "GSM8073222" "GSM8073223" "GSM8073224" "GSM8073225"
# [16] "GSM8073226" "GSM8073227" "GSM8073228" "GSM8073229" "GSM8073230"
# [21] "GSM8073231" "GSM8073232" "GSM8073233" "GSM8073234" "GSM8073235"

all(colnames(expr) %in% rownames(meta))
# [1] FALSE
meta <- meta[colnames(expr), ]
all(colnames(expr) == rownames(meta))
# [1] FALSE

맞습니다. 지금 안 맞는 이유는 양쪽이 서로 다른 샘플 이름 체계를 쓰고 있어서입니다.

expr 열 이름: "4_166__EPG", "20_173__EPG" 같은 내부 샘플 ID
meta 행 이름: "GSM8073211" 같은 GEO accession
meta$title 안에는 "id 166", "id 173" 같은 내부 ID가 들어 있습니다. 또 group은 이미 Control 5 / Other 6 / Responder 14로 잘 만들어져 있습니다.

즉, expr의 166, 173, 145 같은 숫자를 meta$title의 id 166, id 173, id 145와 매칭해야 합니다.

library(stringr)

# 1) expr 열 이름에서 내부 ID 추출: "4_166__EPG" -> "166"
expr_id <- str_match(colnames(expr), "^[0-9]+_([0-9]+)__EPG$")[,2]

# 2) meta$title 에서 내부 ID 추출: "... id 166" -> "166"
meta$id_in_title <- str_match(meta$title, "id\\s+([0-9]+)")[,2]

# 3) expr 순서에 맞게 meta 재정렬
meta2 <- meta[match(expr_id, meta$id_in_title), ]

# 4) 확인
cbind(colnames(expr), expr_id, rownames(meta2), meta2$id_in_title, meta2$title)

#                     expr_id              
#  [1,] "4_166__EPG"  "166"   "NA"    NA NA
#  [2,] "20_173__EPG" "173"   "NA.1"  NA NA
#  [3,] "14_149__EPG" "149"   "NA.2"  NA NA
#  [4,] "18_165__EPG" "165"   "NA.3"  NA NA
#  [5,] "5_168__EPG"  "168"   "NA.4"  NA NA
#  [6,] "16_161__EPG" "161"   "NA.5"  NA NA
#  [7,] "2_146__EPG"  "146"   "NA.6"  NA NA
#  [8,] "9_150__EPG"  "150"   "NA.7"  NA NA
#  [9,] "1_145__EPG"  "145"   "NA.8"  NA NA
# [10,] "21_136__EPG" "136"   "NA.9"  NA NA
# [11,] "13_142__EPG" "142"   "NA.10" NA NA
# [12,] "12_140__EPG" "140"   "NA.11" NA NA
# [13,] "6_138__EPG"  "138"   "NA.12" NA NA
# [14,] "22_137__EPG" "137"   "NA.13" NA NA
# [15,] "8_143__EPG"  "143"   "NA.14" NA NA
# [16,] "15_160__EPG" "160"   "NA.15" NA NA
# [17,] "19_169__EPG" "169"   "NA.16" NA NA
# [18,] "25_158__EPG" "158"   "NA.17" NA NA
# [19,] "7_139__EPG"  "139"   "NA.18" NA NA
# [20,] "17_163__EPG" "163"   "NA.19" NA NA
# [21,] "24_167__EPG" "167"   "NA.20" NA NA
# [22,] "23_155__EPG" "155"   "NA.21" NA NA
# [23,] "11_172__EPG" "172"   "NA.22" NA NA
# [24,] "10_162__EPG" "162"   "NA.23" NA NA
# [25,] "3_152__EPG"  "152"   "NA.24" NA NA


all(expr_id == meta2$id_in_title)
# [1] NA

sum(is.na(meta2$id_in_title))

# [1] 25

좋습니다. 이 결과는 매칭이 아직 안 된 상태라는 뜻입니다.

all(expr_id == meta2$id_in_title)
# [1] NA

sum(is.na(meta2$id_in_title))
# [1] 25

즉,

meta2의 id_in_title이 전부 NA
그래서 match()가 전부 실패
결국 meta$title에서 ID를 제대로 뽑지 못한 것입니다

핵심 원인은 보통 둘 중 하나입니다.

meta$title이 우리가 본 값이 아닌 다른 값으로 들어있다
정규식 추출이 실패했다

먼저 아래를 그대로 실행해 보세요.


colnames(meta)

#  [1] "title"                   "geo_accession"          
#  [3] "status"                  "submission_date"        
#  [5] "last_update_date"        "type"                   
#  [7] "channel_count"           "source_name_ch1"        
#  [9] "organism_ch1"            "characteristics_ch1"    
# [11] "characteristics_ch1.1"   "treatment_protocol_ch1" 
# [13] "growth_protocol_ch1"     "molecule_ch1"           
# [15] "extract_protocol_ch1"    "extract_protocol_ch1.1" 
# [17] "taxid_ch1"               "description"            
# [19] "data_processing"         "data_processing.1"      
# [21] "data_processing.2"       "data_processing.3"      
# [23] "data_processing.4"       "data_processing.5"      
# [25] "data_processing.6"       "platform_id"            
# [27] "contact_name"            "contact_email"          
# [29] "contact_laboratory"      "contact_department"     
# [31] "contact_institute"       "contact_address"        
# [33] "contact_city"            "contact_state"          
# [35] "contact_zip/postal_code" "contact_country"        
# [37] "data_row_count"          "instrument_model"       
# [39] "library_selection"       "library_source"         
# [41] "library_strategy"        "relation"               
# [43] "relation.1"              "supplementary_file_1"   
# [45] "cell line:ch1"           "treatment:ch1"          
# [47] "id_in_title"            


head(meta$title, 10)
 [1] NA NA NA NA NA NA NA NA NA NA

sum(is.na(meta$title))

[1] 25

그리고 ID 추출도 stringr 대신 base R로 더 안전하게 해보겠습니다.

1) expr 쪽 ID 추출

expr_id <- sub("^[0-9]+_([0-9]+)__EPG$", "\\1", colnames(expr))
expr_id
#  [1] "166" "173" "149" "165" "168" "161" "146" "150" "145" "136" "142" "140"
# [13] "138" "137" "143" "160" "169" "158" "139" "163" "167" "155" "172" "162"
# [25] "152"

이건 "4_166__EPG" → "166"으로 바꿉니다.

확인:
head(data.frame(sample = colnames(expr), expr_id = expr_id), 10)

#         sample expr_id
# 1   4_166__EPG     166
# 2  20_173__EPG     173
# 3  14_149__EPG     149
# 4  18_165__EPG     165
# 5   5_168__EPG     168
# 6  16_161__EPG     161
# 7   2_146__EPG     146
# 8   9_150__EPG     150
# 9   1_145__EPG     145
# 10 21_136__EPG     136


2) meta$title 에서 ID 추출
meta$id_in_title <- sub(".*id[[:space:]]+([0-9]+).*", "\\1", meta$title)
head(data.frame(title = meta$title, id_in_title = meta$id_in_title), 10)

#    title id_in_title
# 1   <NA>        <NA>
# 2   <NA>        <NA>
# 3   <NA>        <NA>
# 4   <NA>        <NA>
# 5   <NA>        <NA>
# 6   <NA>        <NA>
# 7   <NA>        <NA>
# 8   <NA>        <NA>
# 9   <NA>        <NA>
# 10  <NA>        <NA>


3) 추출 실패 여부 확인

table(meta$id_in_title == meta$title, useNA = "ifany")

# <NA> 
#   25

만약 대부분 TRUE면, 정규식이 실패해서 원문이 그대로 남은 것입니다.
그 경우 title 형식이 우리가 생각한 것과 다른 겁니다.

그럴 때는 이걸 실행해 보세요.

unique(meta$title)
# [1] NA

4) 매칭 다시 시도
ID가 잘 뽑혔다면:
meta2 <- meta[match(expr_id, meta$id_in_title), ]
rownames(meta2) <- colnames(expr)

all(expr_id == meta2$id_in_title)
# [1] NA
sum(is.na(meta2$id_in_title))
# [1] 25

정상이면:

all(...) → TRUE
sum(is.na(...)) → 0

이 나와야 합니다.


한 번에 실행할 코드

아래를 그대로 복붙해도 됩니다.
# expr 열이름에서 ID 추출
expr_id <- sub("^[0-9]+_([0-9]+)__EPG$", "\\1", colnames(expr))

# meta$title에서 ID 추출
meta$id_in_title <- sub(".*id[[:space:]]+([0-9]+).*", "\\1", meta$title)

# 확인
head(data.frame(sample = colnames(expr), expr_id = expr_id), 10)
head(data.frame(title = meta$title, id_in_title = meta$id_in_title), 10)

# 매칭
meta2 <- meta[match(expr_id, meta$id_in_title), ]
rownames(meta2) <- colnames(expr)

# 점검
all(expr_id == meta2$id_in_title)
sum(is.na(meta2$id_in_title))



만약 여기서도 전부 NA면

그건 meta$title이 실제로는 우리가 본 "tumor, anti-CTLA-4, responder, id 140" 형식이 아니라는 뜻입니다.

그 경우 아래 두 결과만 보내주세요.

head(meta$title, 10)
head(data.frame(sample = colnames(expr)), 10)

# > head(meta$title, 10)
#  [1] NA NA NA NA NA NA NA NA NA NA
# > head(data.frame(sample = colnames(expr)), 10)
#         sample
# 1   4_166__EPG
# 2  20_173__EPG
# 3  14_149__EPG
# 4  18_165__EPG
# 5   5_168__EPG
# 6  16_161__EPG
# 7   2_146__EPG
# 8   9_150__EPG
# 9   1_145__EPG
# 10 21_136__EPG

그럼 제가 바로 딱 맞는 매칭 코드로 바꿔드리겠습니다.

원인 잡혔습니다.

지금 meta$title이 NA로 나온 건, 예전에 했던 이 줄 때문에 meta가 이미 망가졌을 가능성이 큽니다.

meta <- meta[colnames(expr), ]

colnames(expr)는 "4_166__EPG" 같은 이름이고, rownames(meta)는 "GSM8073211" 같은 이름이라 서로 안 맞습니다.
이 상태로 subset하면 내용이 전부 NA인 행들이 만들어집니다. 그 뒤에 meta$title도 NA만 보이게 됩니다.

가장 빠른 해결법

R 세션을 완전히 새로 시작한 뒤, 아래 코드만 순서대로 다시 실행하세요.








range(as.matrix(expr), na.rm = TRUE)

# [1]  Inf -Inf
# Warning messages:
# 1: In min(x, na.rm = na.rm) :
#   no non-missing arguments to min; returning Inf
# 2: In max(x, na.rm = na.rm) :
#   no non-missing arguments to max; returning -Inf

head(expr[,1:5])

#  GSM8073211 GSM8073212 GSM8073213 GSM8073214 GSM8073215





unique(meta$title)

#  [1] "tumor, anti-CTLA-4, responder, id 140"    
#  [2] "tumor, anti-CTLA-4, responder, id 142"    
#  [3] "tumor, anti-CTLA-4, responder, id 149"    
#  [4] "tumor, anti-CTLA-4, responder, id 160"    
#  [5] "tumor, anti-CTLA-4, responder, id 161"    
#  [6] "tumor, anti-CTLA-4, responder, id 163"    
#  [7] "tumor, anti-CTLA-4, responder, id 165"    
#  [8] "tumor, anti-CTLA-4, responder, id 169"    
#  [9] "tumor, anti-CTLA-4, stable, id 173"       
# [10] "tumor, anti-CTLA-4, stable, id 136"       
# [11] "tumor, anti-CTLA-4, stable, id 137"       
# [12] "tumor, anti-CTLA-4, stable, id 155"       
# [13] "tumor, anti-CTLA-4, stable, id 167"       
# [14] "tumor, anti-CTLA-4, stable, id 158"       
# [15] "tumor, anti-CTLA-4, non-responder, id 162"
# [16] "tumor, anti-CTLA-4, non-responder, id 172"
# [17] "tumor, anti-CTLA-4, non-responder, id 138"
# [18] "tumor, anti-CTLA-4, non-responder, id 139"
# [19] "tumor, anti-CTLA-4, non-responder, id 143"
# [20] "tumor, anti-CTLA-4, non-responder, id 150"
# [21] "tumor, IgG2b, control, id 145"            
# [22] "tumor, IgG2b, control, id 146"            
# [23] "tumor, IgG2b, control, id 152"            
# [24] "tumor, IgG2b, control, id 166"            
# [25] "tumor, IgG2b, control, id 168" 


# 👉 여기 결과를 보면 패턴이 보입니다 (예: responder, non-responder 등)