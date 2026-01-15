rm(list = ls())

library(openxlsx)
library(readxl)
# 获取当前工作目录的子文件夹（不包括递归子目录）
sub_dirs <- list.dirs(getwd(), recursive = FALSE)

# 提取所有子文件夹中的.xlsx文件路径
xlsx_files <- list.files(
  path = sub_dirs,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  recursive = FALSE
)

# 将路径向量转换为列表
file_list <- as.list(xlsx_files)

# 打印结果
print(file_list)

# 读取所有文件内容到列表中（每个元素为一个数据框）
data_list <- lapply(file_list, function(file) {
  read_excel(file)  # 读取单个Excel文件
})

# 可选：为列表元素命名（方便后续操作）
names(data_list) <- basename(xlsx_files)

str(data_list)

# 提取所有MR_cell_Results.xlsx的数据
cell_data <- data_list[grepl("MR_cell_Results\\.xlsx", names(data_list))]

# 合并为单个数据框
combined_cell <- do.call(rbind, cell_data)

# 查看合并后的结构
str(combined_cell)

# 提取所有MR_gene_Results.xlsx的数据
gene_data <- data_list[grepl("MR_gene_Results\\.xlsx", names(data_list))]

# 合并为单个数据框
combined_gene <- do.call(rbind, gene_data)

# 查看合并后的结构
str(combined_gene)

phenotype <- read.xlsx("R10_manifest_芬兰  表型下载地址.xlsx",sheet = "Sheet2")

phenotype

# 假设你的数据框名为df
phenotype$phenocode <- paste0("finngen_R10_", phenotype$phenocode)

# 加载包
library(dplyr)

# 合并 cell 数据
combined_cell <- combined_cell %>%
  left_join(phenotype[,1:5], 
            by = c("outcome" = "phenocode"))

# 合并 gene 数据
combined_gene <- combined_gene %>%
  left_join(phenotype[,1:5], 
            by = c("outcome" = "phenocode"))

### ----------------------------------------------------------------------------
###            多重校正
### ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)

# 定义分列校正函数（针对每个变量单独分组校正）
apply_columnwise_correction <- function(data, group_columns) {
  # 遍历每个分组变量
  for (col in group_columns) {
    # 生成校正后的列名
    bh_col <- paste0("pval_BH_", col)
    bonf_col <- paste0("pval_Bonferroni_", col)
    
    # 按当前列分组并校正
    data <- data %>%
      group_by(across(all_of(col))) %>%
      mutate(
        !!bh_col := p.adjust(pval, method = "BH"),
        !!bonf_col := p.adjust(pval, method = "bonferroni")
      ) %>%
      ungroup()
  }
  return(data)
}

# 筛选 nsnp > 3 的行并分档 num_cases
cell_processed <- combined_cell %>%
  filter(nsnp > 3) %>%
  mutate(
    num_cases_group = cut(num_cases,
                          breaks = quantile(num_cases, probs = c(0, 0.33, 0.66, 1)),
                          labels = c("Low", "Medium", "High"),
                          include.lowest = TRUE)
  )

# 定义需要独立校正的列
cell_group_columns <- c("outcome", "exposure", "method", "category", "num_cases_group")

# 应用分列校正
cell_corrected <- apply_columnwise_correction(cell_processed, cell_group_columns)

# 查看新增列
names(cell_corrected)
# 输出示例：
# [1] "pval_BH_outcome"          "pval_Bonferroni_outcome"  "pval_BH_exposure"        
# [4] "pval_Bonferroni_exposure" "pval_BH_method"           ...

gene_processed <- combined_gene %>%
  mutate(
    num_cases_group = cut(num_cases,
                          breaks = quantile(num_cases, probs = c(0, 0.33, 0.66, 1)),
                          labels = c("Low", "Medium", "High"),
                          include.lowest = TRUE)
  )

# 定义需要独立校正的列
gene_group_columns <- c("outcome", "id.exposure", "category", "num_cases_group")

# 应用分列校正
gene_corrected <- apply_columnwise_correction(gene_processed, gene_group_columns)

# 查看新增列
names(gene_corrected)
# 输出示例：
# [1] "pval_BH_outcome"          "pval_Bonferroni_outcome"  "pval_BH_id.exposure"      
# [4] "pval_Bonferroni_id.exposure" "pval_BH_category"        ...

library(openxlsx)

# ==================== 定义通用保存函数 ====================
save_multiple_sheets <- function(data, filename) {
  # 创建新工作簿
  wb <- createWorkbook()
  
  # 添加原始数据分页
  addWorksheet(wb, "Original_Data")
  writeData(wb, sheet = "Original_Data", x = data)
  
  # 获取所有p值列（包括原始pval和校正后的所有p值）
  pval_cols <- grep("^pval", names(data), value = TRUE)
  
  # 遍历每个p值列筛选显著结果
  for (col in pval_cols) {
    # 筛选当前p值列 <0.05 的行
    sig_data <- data[data[[col]] < 0.05, ]
    
    # 仅当存在显著结果时保存
    if (nrow(sig_data) > 0) {
      # 生成分页名称（确保符合Excel规范）
      sheet_name <- substr(col, 1, 31)  # Excel分页名最长31字符
      sheet_name <- gsub("\\.", "_", sheet_name)  # 替换特殊字符
      
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, x = sig_data)
    }
  }
  
  # 保存工作簿
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

# ==================== 执行保存操作 ====================
# 保存 cell 数据
save_multiple_sheets(cell_corrected, "cell_results.xlsx")

# 保存 gene 数据
save_multiple_sheets(gene_corrected, "gene_results.xlsx")

gene_corrected$exposure <- gene_corrected$id.exposure

### ----------------------------------------------------------------------------
###            upset
### ----------------------------------------------------------------------------

# Load required packages
library(openxlsx)
library(UpSetR)
library(dplyr)

# Create output directory
if (!dir.exists("upset result")) {
  dir.create("upset result")
}

