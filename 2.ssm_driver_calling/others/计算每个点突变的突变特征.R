#===============================================================================
#          FILE:Assign_sig_to_single_ssnv.R
# 
#         USAGE:Rscript Assign_sig_to_ssnv.R
# 
#   DESCRIPTION: 这个脚本用于计算986个样本的每个sSNV(jiaoji)的突变特征组成，利用
#                SigProfilerAssignment这个包，可以得到每个突变由各个SBS产生的
#                Probability以及major SBS也就是Probability最大的SBS，最终得到的表格1是每个样本
#                96种突变的major 突变来源sample_id mutation_type Major SBS,表格2是以APOBEC为主
#                要突变来源的突变的SBS Probability分布，Sample，chr，pos，mutation_type,SBS1..，
#                表格3与表格2相同，只是包含所有突变而不只是apobec相关突变
#       OPTIONS: ---
#  REQUIREMENTS: 
#          BUGS: ---
#         NOTES: 
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 11/12/24 
#      REVISION:  ---
#      Reference: DOI: 10.1038/ncomms11383;
#                 https://github.com/AlexandrovLab/SigProfilerAssignmentR
#===============================================================================

#===============================================================================
##1 导入需要的包
if(1){
cat("正在导入包、自写函数和环境变量")
library(tidyverse)
library(reticulate)
library(data.table)
use_python("/home/data/t190513/miniconda3/bin/python")
#py_config()
library(SigProfilerAssignmentR)
#用MutSpot的这个输入就行，也是1-based，start-end改为pos即可
sSNV_filepath <- "/home/data/t190513/1000_noncoding/mutspot/MutSpot/986_jiaoji.mutspot_snv.input"
config_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/wgs_samples_986.txt"
working_dir <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment"
input_dir <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/input"
output_dir <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/output"
result_dir_path <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/Assignment_Solution/Activities"

}
#===============================================================================
##2 导入和预处理文件
if(1){
sSNV <- fread(sSNV_filepath,data.table = F)
names(sSNV) <- c("chr","start","end","ref","alt","sample_id")
config <- read.delim(config_filepath,header = F)
config$V1 <- sapply(config$V1, function(x) {str_remove(x, "(_).*")}) 
config <- as.vector(config$V1) #得到的config是带LC的
if (!file.exists(working_dir)) {
  dir.create(working_dir)
} else {
  cat("目录已存在\n")
}
if (!file.exists(input_dir)) {
  dir.create(input_dir)
} else {
  cat("目录已存在\n")
}
if (!file.exists(output_dir)) {
  dir.create(output_dir)
} else {
  cat("目录已存在\n")
}
}
#===============================================================================
##3 得到目标表格
#3-1 得到输入文件(chr pos(1-based) sample_id ref alt)
for(id in config){
sample_id <- id
input_filename <- paste(sample_id,".vcf",sep ="")
real_input_dir <- file.path(input_dir,"/input")
input_file_path <- file.path(real_input_dir, input_filename)
#如果文件不存在，创建一个
if (!file.exists(input_file_path)) {
df <- sSNV %>% 
  filter(sample_id == {{sample_id}} & chr != "chrM") %>%
  mutate(chr_n=sub("^chr","",chr)) %>%
  select(chr_n,start,sample_id,ref,alt)
write.table(df,file = input_file_path,row.names = F,col.names=F,quote = F,sep = "\t")}}

#3-3 run
exclude_signature_subgroups <- list('UV_signatures','Lymphoid_signatures')
cosmic_fit(samples=input_dir, 
           output=output_dir,
           input_type="vcf",
           context_type="96",
           genome_build="GRCh38",
           cosmic_version=3.3,
           exclude_signature_subgroups=exclude_signature_subgroups,
           export_probabilities_per_mutation=T,
           make_plots = F)

#3-4 得到最终想要的表格
#表格1
#Sample.Names
#MutationType
#Major_SBS
#SBS1
#..
filepath <- file.path(result_dir_path,"Decomposed_MutationType_Probabilities.txt")
df <- read.delim(filepath)
df <- df %>%
  mutate(Major_SBS = apply(df[, 3:73], 1, function(x) {
    # 找出最大值的列索引
    max_col_index <- which.max(x)
    # 获取最大值列的列名
    colnames(df)[max_col_index + 2]  # +2 是因为我们选择了第 3 到第 73 列
  })) %>%
  select(Sample.Names,MutationType,Major_SBS,everything())
output_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/Signature_Single_Mut.txt"
write.table(df,file = output_filepath,row.names = F,quote = F,sep = "\t")
##2 表格2 apobec相关突变
#Sample.Names 样本名
#Chr
#Pos 1-based
#Mutation_Type 96种突变类型的一种
#SBS1 该突变来自SBS1的概率
#SBS2 该突变来自SBS2的概率
#..
#Major_SBS 该突变最有可能来自的SBS
folder_path <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/"
files <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)
data_list <- lapply(files, read.delim)
combined_data <- bind_rows(data_list)
apobec_muts <- combined_data %>%
  filter(SBS2!=0 | SBS13!=0) %>%
  mutate(Major_SBS = apply(combined_data[combined_data$SBS2!=0 | combined_data$SBS13!=0, 5:75], 1, function(x) {
    # 找出最大值的列索引
    max_col_index <- which.max(x)
    # 获取最大值列的列名
    colnames(combined_data)[max_col_index + 4]  # +4 是因为我们选择了第 5 到第 75 列
  })) %>%
  filter(Major_SBS %in% c("SBS2","SBS13"))
apobec_muts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/apobec_muts.txt"
write.table(apobec_muts,file = apobec_muts_filepath,row.names = F,quote = F,sep = "\t")

##表格3
#Sample.Names 样本名
#Chr
#Pos 1-based
#Mutation_Type 96种突变类型的一种
#SBS1 该突变来自SBS1的概率
#SBS2 该突变来自SBS2的概率
#..
#Major_SBS 该突变最有可能来自的SBS

#用到表格中的combined_data变量
all_singlemuts_sigs <- combined_data %>%
  mutate(Major_SBS = apply(combined_data[, 5:75], 1, function(x) {
    # 找出最大值的列索引
    max_col_index <- which.max(x)
    # 获取最大值列的列名
    colnames(combined_data)[max_col_index + 4]  # +4 是因为我们选择了第 5 到第 75 列
  }))
all_singlemuts_sigs_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/all_singlemuts_sigs.txt"
fwrite(all_singlemuts_sigs,file = all_singlemuts_sigs_filepath)
#===============================================================================
#草稿
if(0){
  file_path <- file.path(result_dir_path,"Decomposed_Mutation_Probabilities_1604LC.txt")
  test <- read.
custom_order <- c(1:22, "X", "Y")
df$chr_n <- factor(df$chr_n, levels = custom_order, ordered = TRUE)
df <- df %>%
  arrange(chr_n,start)
files <- list.files(path = result_dir_path, pattern = "*.txt", full.names = TRUE)
data_list <- lapply(files, read.delim)
combined_data <- bind_rows(data_list)
#添加Major_SBS列
cols <- names(combined_data[,5:75])
combined_data <- as.data.table(combined_data)
combined_data[, Major_SBS := apply(.SD, 1, function(row) names(.SD)[which.max(row)], .SDcols = cols )]


data_frame1 <- combined_data %>%
  mutate(Sample_id=Sample.Names,
         Major_SBS=apply(across(SBS1:SBS95), 1, function(row) {
           colnames()[which.max(row)]})) %>%
  select(Sample_id,Chr:SBS95)df <- df %>%
  mutate(Major_SBS = apply(df[, 3:73], 1, function(x) {
    # 查找值为 1 的列的列名
    colnames(df)[which(x == 1) + 2]  # 加 2 因为选择的是第 3 到第 73 列
  })) %>%
  select(Sample.Names,MutationType,Major_SBS,everything())
}