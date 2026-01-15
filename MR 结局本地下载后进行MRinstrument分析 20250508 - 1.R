# 合并重复的库加载（避免重复加载）
library(TwoSampleMR)
library(data.table)
library(tidyverse)
#library(ieugwasr)
library(MRInstruments)
library(mr.raps)
library(MendelianRandomization)
library(vroom)
rm(list = ls())

target_outcome <- c("BRCA") ####根据需求替换

target_gene <- c("CXCL2|CXCL10|CXCL11|CCL5|CCL20|CXCR1|CXCR2|AKT1|MTOR|PIK3CG|PIK3C3") ####根据需求替换,注意必须要1个基因以上才行

#data(gwas_catalog)
#bmi <- subset(gwas_catalog, Phenotype=="Body mass index")
#bmi <- format_data(bmi)

# 设置固定工作目录路径（根据你的需求修改此处）
target_dir <- "E:/20241220 BRCA/20250508/MR"  # 注意使用正斜杠或双反斜杠
setwd(target_dir)  

### 预加载细胞的数据
load("01immune.cell.5e-08.RData")
exposure_dat <- immu.cell.5e.08

# 自动获取ZIP文件并处理路径
list.files(pattern = "\\.gz$")
zip_file <- list.files(pattern = "\\.gz$")[1]          # 取当前目录第一个ZIP文件
stopifnot("未找到gz文件" = !is.na(zip_file))          # 确保文件存在

library(tools)
zip_name <- file_path_sans_ext(zip_file)               # 自动提取文件夹名称

# 0. 清理或创建文件夹（保留文件夹结构）
if (dir.exists(zip_name)) {
  # 清空文件夹内容但保留结构（关键修改）
  unlink(file.path(zip_name, "*"), recursive = TRUE, force = TRUE)
  message("已清空文件夹内容：", zip_name)
} else {
  dir.create(zip_name)
  message("已创建新文件夹：", zip_name)
}

# 1. 创建新文件夹
dir.create(zip_name)                                   # 新建空文件夹

# 结局因素的处理
#读取本地结局
outcome_dat<-vroom(paste0("./",zip_name,".gz"))

# 2. 设置新工作目录
setwd(zip_name)                                        # 切换到解压文件夹
message("当前工作目录已设置为：", getwd())

#head(outcome_dat,10)

library(dplyr)

# 重命名核心列
outcome_dat <- outcome_dat %>%
  rename(
    SNP = rsids,          # SNP标识符
    effect_allele = alt,  # 效应等位基因（如ALT）
    other_allele = ref,   # 非效应等位基因（如REF）
    beta = beta,          # 效应值（原列名已正确）
    se = sebeta,          # 标准误（原列名已正确）
    pval = pval,          # P值（原列名已正确）
    eaf = af_alt          #可选：添加效应等位基因频率列（如果存在）
  )

head(outcome_dat,10)

outcome_dat <- outcome_dat %>% format_data(type = "outcome")
outcome_dat$id.outcome <- zip_name
saveRDS(outcome_dat,"outcome_dat.RDS")

####----------------------------------------------------------------------------
###                基因的分析开始了
####----------------------------------------------------------------------------

### outcome 处理好了
#rm(list = ls())
#outcome_dat <- readRDS("outcome_dat.RDS")
#data("drug_interactions")
data("gtex_eqtl")
summary(as.factor(gtex_eqtl$tissue))

# 筛选CRHBP、MFSD4、MPP7基因
exposure_dat <- subset(gtex_eqtl, 
                       grepl(target_gene, gene_name, ignore.case = TRUE))

summary(as.factor(gtex_eqtl$tissue))

# 在调用format_data前生成唯一ID
exposure_dat$unique_id <- paste(exposure_dat$gene_name, 
                                exposure_dat$tissue, 
                                sep = "_")

head(exposure_dat)

exposure_dat$exposure <- paste(exposure_dat$gene_name,"in",exposure_dat$tissue)
head(exposure_dat)

# 3. 检查每个基因-组织的SNP数量
snp_counts <- exposure_dat %>%
  group_by(gene_name, tissue) %>%
  summarise(N_SNPs = n())

exposure_dat <- format_data(
  exposure_dat,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  id_col = "unique_id",  # 新增的唯一ID列
  eaf_col = "af_alt",     # 假设原始数据包含效应等位基因频率
  gene_col = "gene_name",
  samplesize_col = "n"  # 样本量列（如有）
)

head(exposure_dat)

head(exposure_dat[, c("SNP", "id.exposure", "beta.exposure")])


### harmonise ------------------------------------------------------------------
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
saveRDS(dat,"dat gene.RDS")

library(MRPRESSO)
# 检查有效SNP数量
cat("Available SNPs:", nrow(dat), "\n")
if (nrow(dat) < 10) {
  warning("SNP数量不足，MR-PRESSO无法运行！建议放宽筛选条件。")
} else {
  # 运行MR-PRESSO的代码
  # MR-PRESSO全局检验（需安装MRPRESSO包）
  presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                BetaExposure = "beta.exposure",
                                SdOutcome = "se.outcome",
                                SdExposure = "se.exposure",
                                data = dat,
                                SignifThreshold = 0.05,  # 默认0.05，放宽至0.1
                                NbDistribution = 1000)# 默认1000
  
  cat("MR-PRESSO global test p-value:", presso$`Global Test`$Pvalue, "\n")
}

MR_PRESSO <- as.data.frame(presso)

### MR 分析--------------------------------------------------------------------
mr_method_list <- mr_method_list()

print(mr_method_list$obj)

# MR 分析
result <- mr(dat, method_list = mr_method_list[-16,]$obj)#去掉mr raps

# 计算OR值
result <- generate_odds_ratios(result)

result$outcome <- zip_name

head(result)

# 异质性检验
mr_heterogeneity(dat)
#若Q检验P < 0.05：存在显著异质性，提示潜在水平多效性。

# 水平多效性检验
mr_pleiotropy_test(dat)

# 散点图
p1 <- mr_scatter_plot(result, dat)

# 森林图
result_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(result_single)

# 留一图
result_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(result_loo)

# 漏斗图
result_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(result_single)

# 安装并加载所需包（如果尚未安装）
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# 创建 Excel 工作簿对象
wb <- createWorkbook()

# 添加第一个工作表：主要结果（result）
addWorksheet(wb, sheetName = "Main Results")
writeData(wb, sheet = 1, x = result, startCol = 1, startRow = 1)

# 添加第二个工作表：单SNP分析（result_single）
addWorksheet(wb, sheetName = "Single SNP")
writeData(wb, sheet = 2, x = result_single, startCol = 1, startRow = 1)

# 添加第三个工作表：留一法分析（result_loo）
addWorksheet(wb, sheetName = "Leave-One-Out")
writeData(wb, sheet = 3, x = result_loo, startCol = 1, startRow = 1)

# 添加第4个工作表：MRPRESSO
addWorksheet(wb, sheetName = "MRPRESSO")
writeData(wb, sheet = 4, x = MR_PRESSO, startCol = 1, startRow = 1)

# 保存Excel文件到当前工作目录
saveWorkbook(wb, 
             file = "MR_gene_Results.xlsx",
             overwrite = TRUE)

# 确认输出完成
message("Results saved to: ", file.path(getwd(), "MR_gene_Results.xlsx"))


####----------------------------------------------------------------------------
###                细胞的分析开始了
####----------------------------------------------------------------------------
exposure_dat <- immu.cell.5e.08

### harmonise ------------------------------------------------------------------
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,
                      action = 3)  # action=3时自动翻转冲突方向
saveRDS(dat,"dat cell.RDS")

mr_method_list <- mr_method_list()

print(mr_method_list$obj)

#    TwoSampleMR   包默认使用
#五种常用方法（IVW、MR-Egger、Weighted Median、Simple Mode、Weighted Mode）
#进行孟德尔随机化分析
# MR 分析
result <- mr(dat
             #, method_list = mr_method_list[-16,]$obj
)#去掉mr raps

# 计算OR值
result <- generate_odds_ratios(result)

result$outcome <- zip_name

head(result)

# 异质性检验
mr_heterogeneity(dat)

# 水平多效性检验
mr_pleiotropy_test(dat)

# 散点图
p1 <- mr_scatter_plot(result, dat)

# 森林图
result_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(result_single)

# 留一图
result_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(result_loo)

# 漏斗图
result_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(result_single)



head(result)
head(result_loo)
head(result_single)

# 通过GitHub安装最新开发版（修复潜在bug）
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("rondolab/MR-PRESSO",force = TRUE)

# 安装并加载所需包（如果尚未安装）
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# 创建 Excel 工作簿对象
wb <- createWorkbook()

# 添加第一个工作表：主要结果（result）
addWorksheet(wb, sheetName = "Main Results")
writeData(wb, sheet = 1, x = result, startCol = 1, startRow = 1)

# 添加第二个工作表：单SNP分析（result_single）
addWorksheet(wb, sheetName = "Single SNP")
writeData(wb, sheet = 2, x = result_single, startCol = 1, startRow = 1)

# 添加第三个工作表：留一法分析（result_loo）
addWorksheet(wb, sheetName = "Leave-One-Out")
writeData(wb, sheet = 3, x = result_loo, startCol = 1, startRow = 1)

# 保存Excel文件到当前工作目录
saveWorkbook(wb, 
             file = "MR_cell_Results.xlsx",
             overwrite = TRUE)

# 确认输出完成
message("Results saved to: ", file.path(getwd(), "MR_cell_Results.xlsx"))
