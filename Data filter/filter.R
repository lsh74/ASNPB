rm(list = ls())
library(dplyr)
library(tidyverse)
setwd("/data3/lishuhan/GSEdir/")
combined_data<-read.csv("./filter_data.csv")
############年龄信息筛选################################################################

# 提取列名
col_names <- names(combined_data)

# 提取包含关键词的列名及其对应的列号
age_cols <- grep("*(age|Age)", col_names, value = TRUE)
age_indices <- grep("*(age|Age)", col_names)
age_info <- paste(age_indices, age_cols, sep = ": ")

# 将列名及其对应的列号写入txt文件
writeLines(age_info, "./result/age_cols.txt")

####检查不确定的列
therapy_non <- subset(combined, !is.na(X.Sample_characteristics_ch1..biological.sex.))

# 进一步检查并过滤空字符串或仅包含空格的字符串
therapy_values <- subset(therapy_non, nzchar(trimws(X.Sample_characteristics_ch1..biological.sex.)))
therapy_values$X.Sample_characteristics_ch1..biological.sex.



##########将列名里指出年龄是以月份、天来计数的列提出来修改以年为单位,并放回combined_data中##########
column_ind <- c(539,925,1326,788,1625,1956,2096)
extracted_columns <- combined_data[, column_ind]
has_non_empty_value <- apply(extracted_columns, 1, function(row) {
  any(row != "" & !is.na(row))
})
# 保留至少有一个非空值的行
filtered <- extracted_columns[has_non_empty_value,]

#####仅保留至少包含以上某一列内容的行
# 检查每一列是否包含非数字内容
for (i in 1:ncol(filtered)) {
  column <- filtered[, i]
  
  # 检查列中的每个元素是否是数字
  is_numeric <- sapply(column, function(x) {
    # 如果元素是空的或者是 NA，不考虑它
    if (x == "" || is.na(x)) {
      return(TRUE)
    }
    # 检查元素是否是数字
    return(is.numeric(as.numeric(x)) && !is.na(as.numeric(x)))
  })
  
  # 如果有任何非数字内容，打印列索引
  if (!all(is_numeric)) {
    cat("Column index:", i, "\n")
  }
}
filtered#1539
write.csv(filtered,"./result/m_d_age_filtered.csv")

# 读取数据
file_path <- "./result/m_d_age_filtered.csv"
data <- read.csv(file_path, row.names = 1)

# 转换数据到年为单位
data_in_years <- data
data_in_years[,-4] <- data_in_years[,-4] / 12  # 对于非第五列（排除行名），除以12
data_in_years[, 4] <- data_in_years[, 4] / 365 # 对于第五列（排除行名），除以365

# 四舍五入到两位小数
data_in_years <- round(data_in_years, 2)

# 保存转换后的数据
write.csv(data_in_years, "./result/m_age_filtered.csv", row.names = TRUE)

###合回去combined——data
# 读取 modified_data 和 combined_data
modified_file_path <- "./result/m_age_filtered.csv"
modified_data <- read.csv(modified_file_path, row.names = 1)

# 获取行名和列名
row_names <- rownames(modified_data)
col_names <- colnames(modified_data)

# 用 modified_data 中的数据替换 combined_data 中对应的数据
for (row in row_names) {
  for (col in col_names) {
    if (row %in% rownames(combined_data) && col %in% colnames(combined_data)) {
      combined_data[row, col] <- modified_data[row, col]
    }
  }
}

# 保存更新后的 combined_data 到一个新的CSV文件
write.csv(combined_data, "./result/updated_combined_data.csv", row.names = TRUE)
###################################将所有包含年龄信息的列提出来，修改之后再放入combined_data中####################
rm(list = ls())
combined<-read.csv("./age_reslut/updated_combined_data.csv")
##################对每列的结果分别存入不同csv文件，第一列是非数字的内容，第二列是对应在filtered_rows的行号也是conbind的行号
cols_to_keep <- c(15, 18, 85, 202, 214, 222, 232, 290, 357, 376, 382, 431, 432, 450, 453, 539, 
                  561, 608, 640, 771, 777, 781, 788, 815, 820, 825, 832, 866, 892, 925, 961, 962, 
                  985, 986, 988, 1178, 1201, 1205, 1215, 1326, 1385, 1430, 1455, 1461, 1505, 1579,
                  1594, 1625, 1660, 1663, 1763, 1809, 1898, 1956, 1958, 2031, 2057, 2078, 2096, 
                  2219, 2295, 2337, 2414, 2443)
filtered_rows <- combined[, cols_to_keep]
# 第一部分：创建包含非数字内容的CSV文件

# 获取 filtered_rows 的行名，以便找到在 combined_data 中的原始行号
filtered_row_names <- rownames(filtered_rows)

for (i in 1:ncol(filtered_rows)) {
  column <- filtered_rows[, i]
  
  # 获取在 combined_data 中的原始列索引
  original_col_index <- cols_to_keep[i]
  
  # 获取当前列的名字
  column_name <- colnames(combined)[original_col_index]
  
  # 找到非数字内容的索引
  non_numeric_indices <- which(sapply(column, function(x) {
    # 如果元素是空的或者是 NA，不考虑它
    if (x == "" || is.na(x)) {
      return(FALSE)
    }
    # 判断元素是否不是数字
    return(!is.numeric(as.numeric(x)) || is.na(as.numeric(x)))
  }))
  
  # 获取在 combined_data 中的原始行号
  original_row_numbers <- filtered_row_names[non_numeric_indices]
  
  # 如果存在非数字内容
  if (length(non_numeric_indices) > 0) {
    # 创建一个数据框，包含非数字内容及其原始行号
    non_numeric_contents <- data.frame(
      'Content' = column[non_numeric_indices],
      'Row_Number' = original_row_numbers
    )
    
    # 使用 combined_data 的列名创建文件名
    file_name <- paste0("./GSEdir/age_reslut/", column_name, ".csv")
    
    # 将数据框写入 CSV 文件
    write.csv(non_numeric_contents, file_name, row.names = FALSE)
  }
}


# 第二部分：加载CSV文件并替换 combined_data 中的内容

# 获取CSV文件的路径
csv_files <- list.files("./age_reslut/cols/", full.names = TRUE)

# 遍历每个CSV文件
for (file_path in csv_files) {
  # 读取CSV文件
  csv_data <- read.csv(file_path)
  
  # 获取列名（从文件名解析）
  column_name <- gsub(".*\\/|\\.csv", "", file_path)
  
  # 遍历每行数据并替换相应的内容
  for (i in 1:nrow(csv_data)) {
    content <- csv_data$Content[i]
    row_number <- as.numeric(csv_data$Row_Number[i])
    
    # 替换 combined_data 中的数据
    combined[row_number, column_name] <- content
  }
}

# 保存修改后的 combined_data
write.csv(combined, "./age_reslut/linshi/change_age_col_combined.csv", row.names = FALSE)
################将年龄的列并成一列age###################
combined<-read.csv("./age_reslut/linshi/change_age_col_combined.csv")
columns_to_merge <- c(15, 18, 85, 202, 214, 222, 232, 290, 357, 376, 382, 431, 432, 450, 453, 539,
                      561, 608, 640, 771, 777, 781, 788, 815, 820, 825, 832, 866, 892, 925, 961, 962,
                      985, 986, 988, 1178, 1201, 1205, 1215, 1326, 1385, 1430, 1455, 1461, 1505, 1579,
                      1594, 1625, 1660, 1663, 1763, 1809, 1898, 1956, 1958, 2031, 2057, 2078, 2096, 
                      2219, 2295, 2337, 2414, 2443)
# 为所有这些列中的空字符串（""）赋值为 NA
combined[columns_to_merge] <- lapply(combined[columns_to_merge], function(x) ifelse(x=="", NA, x))

# 将所有要合并的列转换为字符类型
combined[columns_to_merge] <- lapply(combined[columns_to_merge], as.character)

# 合并指定的列
combined <- combined %>%
  mutate(age = dplyr::coalesce(!!!.[columns_to_merge]))

# 删除合并后的原始列
combined <- combined[, setdiff(names(combined), names(combined[columns_to_merge]))]

# 删除 age 列中全是 NA 或空白的行
combined <- combined %>%
  dplyr::filter(!is.na(age) & age != "")

# 检查并删除所有全是 NA 或空白的列
empty_columns <- sapply(combined, function(x) all(is.na(x) | x == ""))
combined <- combined[, !empty_columns]
combined <- dplyr::select(combined, 1:2, age, everything())

# 写入到新的 CSV 文件
write.csv(combined, "./age_reslut/linshi/up_age_data.csv", row.names = FALSE)#82473*1465

#############筛选性别信息---------------------------------------------------------------------------------
col_names <- names(combined)
# 包含"gender"或"sex"的列
gender_cols <- grep("*(gender|sex|Gender|Sex|SEX|GENDER|male|female|gendr)", col_names, value = TRUE)
gender_indices <- grep("*(gender|sex|Gender|Sex|SEX|GENDER|male|female|gendr)", col_names)
gender_info <- paste(gender_indices, gender_cols, sep = ": ")
writeLines(gender_info, "./result/sex_cols.txt")
#####年龄信息&性别信息（可选）------------
combined <- read.csv("./age_reslut/linshi/up_age_data.csv")

# 指定列进行替换
cols_0_male_1_female <- c(96,543,940,1176)
cols_0_female_1_male <- c(517,865)
col_1_male_2_female <- 838

# 执行替换
combined[cols_0_male_1_female] <- lapply(combined[cols_0_male_1_female], function(x) ifelse(x == 0, 'male', ifelse(x == 1, 'female', x)))
combined[cols_0_female_1_male] <- lapply(combined[cols_0_female_1_male], function(x) ifelse(x == 0, 'female', ifelse(x == 1, 'male', x)))
combined[col_1_male_2_female] <- ifelse(combined[[col_1_male_2_female]] == 1, 'male', 
                                        ifelse(combined[[col_1_male_2_female]] == 2, 'female', combined[[col_1_male_2_female]]))


# 指定列索引用于合并
columns_to_merge <- c(14, 75, 83, 96, 153, 281, 285, 320, 322, 456, 517, 543,709,838,865,881,940,1165,1176,1249,1312,1399,1462)

# 为所有这些列中的空字符串（""）赋值为 NA
combined[columns_to_merge] <- lapply(combined[columns_to_merge], function(x) ifelse(x=="", NA, x))

# Convert all columns to character before combining
combined[columns_to_merge] <- lapply(combined[columns_to_merge], as.character)

# Merge columns into single 'sex' column
combined <- combined %>%
  rowwise() %>%
  mutate(sex = paste(na.omit(c_across(all_of(columns_to_merge))), collapse = "")) %>%
  ungroup() %>%
  dplyr::select(-all_of(columns_to_merge))

# 将sex列提到第4列
column_order <- c(1:3, ncol(combined), 4:(ncol(combined)-1))
combined <- combined[, column_order]
write.csv(combined, "./sex_result/2up_age_sex_data.csv", row.names = FALSE)#82473 * 1443


###RNA_seq筛选---------------------
combined <- combined[combined$X.Sample_library_strategy=='RNA-Seq',]#23127 * 1443
write.csv(combined, "./seq_reslut/2up_age_sex_rna_data.csv", row.names = FALSE)#23127 * 1443
#######逐个检查GSE---------------------------------------------------------------------------
combined<-read.csv("./seq_reslut/2up_age_sex_rna_data.csv")
# 提取第50列的唯一值

library(readr)

# 读取 CSV 文件
data <- read.csv("./result/final.csv")

# 指定你的 GSE 号列表
gse_list <- c("GSE232186", "GSE230370", "GSE218212", "GSE228702", "GSE218843", "GSE218213", "GSE228320", "GSE222393", "GSE227743", "GSE224714", "GSE225217",
              "GSE223885", "GSE213912", "GSE220969", "GSE199750", "GSE212584", "GSE205402", "GSE221911", "GSE142760", "GSE221091", "GSE200462", "GSE210943",
              "GSE220589", "GSE216030", "GSE141529", "GSE205867", "GSE214137", "GSE213004", "GSE214240", "GSE214207", "GSE213486", "GSE185011", "GSE212038",
              "GSE197903", "GSE207057", "GSE209591", "GSE176216", "GSE195796", "GSE196117", "GSE153315", "GSE202625", "GSE201754", "GSE199230", "GSE198856",
              "GSE198123", "GSE204864", "GSE198449", "GSE193952", "GSE181179", "GSE165080", "GSE195599", "GSE189126", "GSE189125", "GSE186144", "GSE169381",
              "GSE182031", "GSE178408", "GSE185263", "GSE185863", "GSE190510", "GSE180387", "GSE178967", "GSE173697", "GSE174083", "GSE186928", "GSE185058",
              "GSE142444", "GSE169687", "GSE154616", "GSE176153", "GSE175701", "GSE152693", "GSE175604", "GSE166663", "GSE173432", "GSE155897", "GSE157967",
              "GSE148885", "GSE157344", "GSE161731", "GSE163219", "GSE152641", "GSE161031", "GSE133397", "GSE134637", "GSE138544", "GSE136849", "GSE157148",
              "GSE149729", "GSE159121", "GSE151161", "GSE123658", "GSE142844", "GSE148986", "GSE125873", "GSE130279", "GSE143507", "GSE143567", "GSE135936",
              "GSE139242", "GSE134080", "GSE102114", "GSE116672", "GSE122485", "GSE119193", "GSE117769", "GSE116899", "GSE94438",  "GSE112087", "GSE104423",
              "GSE103232", "GSE110041", "GSE75337",  "GSE107196", "GSE85263",  "GSE107981", "GSE79362",  "GSE51799",  "GSE52166")

# 提取第二列 GSE 号在列表中的行
filtered_data <- data[!data[[2]] %in% gse_list,]

# 打印筛选出的数据
print(filtered_data)
write.csv(filtered_data,"./result/addfinal.csv")
# 打印唯一值
print(unique_values)
write.csv(unique_values,"./result/final.csv")
# [1] "GSE232186" "GSE230370" "GSE218212" "GSE228702" "GSE218843" "GSE218213" "GSE228320" "GSE222393" "GSE227743" "GSE224714" "GSE225217"
# [12] "GSE223885" "GSE213912" "GSE220969" "GSE199750" "GSE212584" "GSE205402" "GSE221911" "GSE142760" "GSE221091" "GSE200462" "GSE210943"
# [23] "GSE220589" "GSE216030" "GSE141529" "GSE205867" "GSE214137" "GSE213004" "GSE214240" "GSE214207" "GSE213486" "GSE185011" "GSE212038"
# [34] "GSE197903" "GSE207057" "GSE209591" "GSE176216" "GSE195796" "GSE196117" "GSE153315" "GSE202625" "GSE201754" "GSE199230" "GSE198856"
# [45] "GSE198123" "GSE204864" "GSE198449" "GSE193952" "GSE181179" "GSE165080" "GSE195599" "GSE189126" "GSE189125" "GSE186144" "GSE169381"
# [56] "GSE182031" "GSE178408" "GSE185263" "GSE185863" "GSE190510" "GSE180387" "GSE178967" "GSE173697" "GSE174083" "GSE186928" "GSE185058"
# [67] "GSE142444" "GSE169687" "GSE154616" "GSE176153" "GSE175701" "GSE152693" "GSE175604" "GSE166663" "GSE173432" "GSE155897" "GSE157967"
# [78] "GSE148885" "GSE157344" "GSE161731" "GSE163219" "GSE152641" "GSE161031" "GSE133397" "GSE134637" "GSE138544" "GSE136849" "GSE157148"
# [89] "GSE149729" "GSE159121" "GSE151161" "GSE123658" "GSE142844" "GSE148986" "GSE125873" "GSE130279" "GSE143507" "GSE143567" "GSE135936"
# [100] "GSE139242" "GSE134080" "GSE102114" "GSE116672" "GSE122485" "GSE119193" "GSE117769" "GSE116899" "GSE94438"  "GSE112087" "GSE104423"
# [111] "GSE103232" "GSE110041" "GSE75337"  "GSE107196" "GSE85263"  "GSE107981" "GSE79362"  "GSE51799"  "GSE52166"

################所有可用健康样本整合，画图展示---------------------------------
#GSE228702,GSE197903,GSE209591,GSE196117,GSE153315
# GSE202625,GSE198123,GSE193952,GSE195599,GSE173697
# GSE169687,GSE176153,GSE161731,GSE163219,GSE134637
# GSE123658,GSE134080,GSE102114,GSE94438,GSE112087
# GSE110041,GSE79362,GSE51799,GSE124326,GSE181228
# GSE136371,GSE122485,GSE222889,GSE191238,GSE182038,GSE120312

################################所有能用的样本信息统计--------------------------------
library(dplyr)
setwd("./")
data<-read.csv("/data3/lishuhan/GSEdir/final_data/GSMinfo/upcombined.csv")
# [1] "GSE51799"  "GSE79362"  "GSE94438"  "GSE102114" "GSE110041" "GSE112087" "GSE122485" "GSE123658" "GSE124326" "GSE134080" "GSE134637" "GSE136371" "GSE152641"
# [14] "GSE153315" "GSE161731" "GSE166663" "GSE169687" "GSE176153" "GSE181228" "GSE193952" "GSE195599" "GSE196117" "GSE197903" "GSE198123" "GSE202625" "GSE209591"
# [27] "GSE223885" "GSE120312" "GSE182038" "GSE191238" "GSE222889" "GSE173697" "GSE228702"
###合并年龄------------------------------------------------------------------
columns_to_merge <- c(13,124,134,137,139,141)


data[columns_to_merge] <- lapply(data[columns_to_merge], function(x) ifelse(x=="", NA, x))


# 将所有要合并的列转换为字符类型
data[columns_to_merge] <- lapply(data[columns_to_merge], as.character)

# 合并指定的列
data <- data %>%
  mutate(age = dplyr::coalesce(!!!.[columns_to_merge]))

# 获取最后一列的名称
last_column_name <- names(data)[ncol(data)]

# 移动最后一列到第三列
data <- data %>%
  dplyr::select(1:2, last_column_name, everything())

library(dplyr)
library(stringr)

# 移除 " years"，只保留数字
data <- data %>%
  mutate(age = str_replace(age, " years", ""))

# 删除包含"NA"，"post 20 years old"，"not collected"的行
data <- data %>%
  dplyr::filter(!age %in% c("NA", "post 20 years old", "not collected"))

data <- data %>%
  dplyr::slice(-c(917))

data <- data %>%
  dplyr::slice(-c(1354:1360))

data$age[c(1420,1422)] <- 73





write.csv(data,"./final_data/GSMinfo/2up_age_data.csv")1550*1447

###合并性别---------------------------------------------------------------------
cols_0_male_1_female <- c(15)
cols_0_female_1_male <- c(112)
data[cols_0_male_1_female] <- lapply(data[cols_0_male_1_female], function(x) ifelse(x == 0, 'male', ifelse(x == 1, 'female', x)))
data[cols_0_female_1_male] <- lapply(data[cols_0_female_1_male], function(x) ifelse(x == 0, 'female', ifelse(x == 1, 'male', x)))


col_names <- names(data)
ender_cols <- grep("*(gender|sex|Gender|Sex|SEX|GENDER|male|female|gendr)", col_names, value = TRUE)
gender_indices <- grep("*(gender|sex|Gender|Sex|SEX|GENDER|male|female|gendr)", col_names)
gender_info <- paste(gender_indices, ender_cols, sep = ": ")
writeLines(gender_info, "./result/gender_sex_cols.txt")

# 寻找包含特定字符串的列的索引
columns_to_merge <- grep("gender|sex|Gender|Sex|SEX|GENDER|male|female|gendr", names(data))

# 转换这些列为字符类型
data[columns_to_merge] <- lapply(data[columns_to_merge], as.character)

# 合并这些列
gender_cols <- grep("gender|sex|Gender|Sex|SEX|GENDER|male|female|gendr", names(data), value = TRUE)
data$sex <- do.call(dplyr::coalesce, data[gender_cols])


cols_to_transform <- c("X.Sample_characteristics_ch1..gender..0.male..1.female..",
                       "X.Sample_characteristics_ch1..gender.",
                       "X.Sample_characteristics_ch1..Sex.",
                       "X.Sample_characteristics_ch1..sex..0.females..1.male..")

# 为所有这些列中的空字符串（""）赋值为 NA
data[cols_to_transform] <- lapply(data[cols_to_transform], function(x) ifelse(x=="", NA, x))

# 合并这些列为 "sex" 列
data$sex <- do.call(dplyr::coalesce, data[cols_to_transform])

# 替换 "sex" 列中的 "M" 为 "male"，"F" 为 "female", "Male"为"male", "Female"为"female"
data$sex <- ifelse(data$sex %in% c("M","MALE","Male"), "male",
                   ifelse(data$sex %in% c("F", "FEMALE","Female"), "female", data$sex))

# 将 "sex" 列移到第四列
last_column_name <- names(data)[ncol(data)]

data <- data %>%
  dplyr::select(1:3, last_column_name, everything())

write.csv(data,"./final_data/GSMinfo/uphealth_data.csv")#1543*148



values_to_remove <- c("GSE152641")
rows_to_keep <- !(data[, 41] %in% values_to_remove)
data <- data[rows_to_keep, ]
write.csv(data,"./final_data/GSMinfo/2uphealth_data.csv")#1519*148


#####最终的最终
# "GSE51799"  "GSE79362"  "GSE94438"  "GSE102114" "GSE110041" 
# "GSE112087" "GSE122485" "GSE123658" "GSE124326" "GSE134080" 
# "GSE134637" "GSE136371" "GSE153315""GSE161731" "GSE169687" 
# "GSE176153" "GSE181228" "GSE193952" "GSE195599" "GSE196117" 
# "GSE197903" "GSE198123" "GSE202625" "GSE209591" "GSE120312" 
# "GSE182038" "GSE191238" "GSE222889" "GSE173697" "GSE228702"
###整理一下1255个样本信息---------------------------
###GSE94438
df2_subset <- a[, c(1,12,18,38)]

# 合并数据框，使用 Sample_Name 和 X.SAMPLE 作为连接键
merged_df <- merge(df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
colnames(merged_df)[8:10] <- c("Tissue", "Phenotype", "Data.Type")

##GSE79362
G7<-read.csv("/data3/lishuhan/GSEdir/final30/GSE79362/combined.csv")
df2_subset <- G7[, c(1,20,12,37)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )

# 删除多余的列
merged_df<-merged_df[,-(8:13)]

##GSE102114
G1<-read.csv("/data3/lishuhan/GSEdir/final30/GSE102114/combined.csv")
df2_subset <- G1[, c(1,13,12,34)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )

# 删除多余的列
merged_df<-merged_df[,-(8:13)]

##GSE122485
G12<-read.csv("/data3/lishuhan/GSEdir/final30/GSE122485/combined.csv")
df2_subset <- G12[, c(1,12,15,35)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )

# 删除多余的列
merged_df<-merged_df[,-(8:13)]

##GSE123658
G2<-read.csv("/data3/lishuhan/GSEdir/final30/GSE123658/combined.csv")
df2_subset <- G2[, c(1,12,13,36)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE124326
G3<-read.csv("/data3/lishuhan/GSEdir/final30/GSE124326/combined.csv")
df2_subset <- G3[, c(1,18,12,64)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE134080
G4<-read.csv("/data3/lishuhan/GSEdir/final30/GSE134080/combined.csv")
df2_subset <- G4[, c(1,16,12,35)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE134637
G5<-read.csv("/data3/lishuhan/GSEdir/final30/GSE134637/combined.csv")
df2_subset <- G5[, c(1,12,17,35)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
df2_subset$Phenotype<-"Healthy"
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE136371
G6<-read.csv("/data3/lishuhan/GSEdir/final30/GSE136371/combined.csv")
df2_subset <- G6[, c(1,13,12,34)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE161731
G8<-read.csv("/data3/lishuhan/GSEdir/final30/GSE161731/combined.csv")
df2_subset <- G8[, c(1,9,16,37)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE169687
G9<-read.csv("/data3/lishuhan/GSEdir/final30/GSE169687/combined.csv")
df2_subset <- G9[, c(1,9,19,36)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE173697
G10<-read.csv("/data3/lishuhan/GSEdir/final30/GSE173697/combined.csv")
df2_subset <- G10[, c(1,9,12,32)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE181228
G11<-read.csv("/data3/lishuhan/GSEdir/final30/GSE181228/combined.csv")
df2_subset <- G11[, c(1,9,12,32)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE182038
G13<-read.csv("/data3/lishuhan/GSEdir/final30/GSE182038/combined.csv")
df2_subset <- G13[, c(1,9,15,34)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
df2_subset$Phenotype<-"Healthy"
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE191238
G14<-read.csv("/data3/lishuhan/GSEdir/final30/GSE191238/combined.csv")
df2_subset <- G14[, c(1,9,13,35)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE193952
G15<-read.csv("/data3/lishuhan/GSEdir/final30/GSE193952/combined.csv")
df2_subset <- G15[, c(1,13,2,33)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE195599
G16<-read.csv("/data3/lishuhan/GSEdir/final30/GSE195599/combined.csv")
df2_subset <- G16[, c(1,9,12,34)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")

merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE197903
G18<-read.csv("/data3/lishuhan/GSEdir/final30/GSE197903/combined.csv")
df2_subset <- G18[, c(26,28,14,3)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
colnames(df2_subset)[1]<-"X.SAMPLE"
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE198123
G19<-read.csv("/data3/lishuhan/GSEdir/final30/GSE198123/combined.csv")
df2_subset <- G19[, c(1,9,13,31)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE209591
G20<-read.csv("/data3/lishuhan/GSEdir/final30/GSE209591/combined.csv")
df2_subset <- G20[, c(1,9,2,33)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]

##GSE222889
G21<-read.csv("/data3/lishuhan/GSEdir/final30/GSE222889/combined.csv")
df2_subset <- G21[, c(1,9,13,33)]
colnames(df2_subset)[2:4] <- c("Tissue", "Phenotype", "Data.Type")
merged_df <- merge(merged_df, df2_subset, by.x = "Sample_Name", by.y = "X.SAMPLE", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    Phenotype = coalesce(Phenotype.x, Phenotype.y),
    Data.Type = coalesce(Data.Type.x, Data.Type.y)
  )
merged_df<-merged_df[,-(8:13)]
write.csv(merged_df,"/data3/lishuhan/Article/附表/sample_info.csv")