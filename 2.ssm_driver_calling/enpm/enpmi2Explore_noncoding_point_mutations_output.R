#===============================================================================
#
#          FILE:enpmi2Explore_noncoding_point_mutations_output.R
# 
#         USAGE: ./enpmi2Explore_noncoding_point_mutations_output.R
# 
#   DESCRIPTION: 这个脚本用于分析986fdscc突变频率超过10的点突变及其与靶基因表达的关系
#       OPTIONS: ---
#  REQUIREMENTS: 
#          BUGS: ---
#         NOTES: 目前使用的输入是986 jiaoji hg38
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 1/16/25
#      REVISION:  ---
#      Reference:
#===============================================================================

#输入是vep注释后的文件
#0 导入需要的包和自写函数和路径,设置环境变量
if(1)
{
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(data.table)
  library(scales)
  library(ggbreak)
  library(stringr)
  library(GenomicRanges)
  library(survminer)
  library(survival)
  library(doParallel)
  library(ComplexHeatmap)
  library(fdrtool)
  library(RColorBrewer)
config_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/wgs_samples_986.txt"
cancer_gene_737_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/737 cancer genes by pcawg.csv"
lung_drivers_df_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/20220203_lung_drivers.csv"
cosmic_census_cancer_gene_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/Cancer Gene Census v100.tsv"
enpmi_filepath <- "/home/data/t190513/1000_noncoding/nc_point_mutations/986_jiaoji.enpm.input"
annotations_hg38_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_elements_hg38_noctcf.bed"
TAD_withgenes_hg38_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/A549_TAD_withgenes_hg38.bed"
normalized_counts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/counts_1000_log2tpmplus1.txt"
clinical_data_file_path <- "/home/data/t190513/1000_noncoding/activedriverwgs/clin_info_1000_with_pathway_with_gene.txt"
all_gene_bed_hg38_file_path <- "/home/data/t190513/1000_noncoding/nc_point_mutations/all_gene_hg38.bed4"
oncokb_cancergene_list_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/oncokb_cancerGeneList_2024_10_24.tsv"
cnv_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_data_by_genes.txt"
GEM_mappability_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/GRCh38_mappability_100mer.gem.bed"
Encode_blacklist_regions_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/encode_hg38.blacklist.bed"
IR_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/IR_hg38.rds"
all_singlemuts_sigs_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/all_singlemuts_sigs.txt"
folder_path <- "/home/data/t190513/1000_noncoding/nc_point_mutations"
#Func-1 输入一个或多个以;分隔的基因名(Ensemble_id或Hugo_Symbol都行，混着也行),其突变的样本集合(Tumor_Sample_Barcode,逗号相连的一串值),归一化矩阵的行名以及Hugo_symbol
#返回突变组相对非突变组的表达p值
p_expr <- function(genes,mutated_samples,nc_rownames,nc_hugosymbol){
  #如果genes不是NA，开弄
  if(all(!is.na(genes))){
  genes <- strsplit(genes,";")[[1]]
  mutated_samples <- strsplit(mutated_samples,",")[[1]]
  result <- rep(NA,length(genes))
  # 预先获取mutated_samples对应的列索引
  mutated_indices <- which(names(normalized_counts) %in% mutated_samples)
  for(i in 1:length(genes)){
  #首先匹配到归一化counts矩阵该gene对应的行
  d_f <- data.frame()
  #如果该基因是ENSG开头,与在rownames里即可
  if(grepl("^ENSG",genes[i])){
    d_f <- normalized_counts[grepl(genes[i],nc_rownames),]
  }else{ #如果是Hugo_Symbol与Hugo_Symbol比对，需完全匹配
    d_f <- normalized_counts[nc_hugosymbol==genes[i],]
  }
  #如果能匹配上
  if(nrow(d_f)==1){
  # 将表达量分组
  #原本的mutated samples是以，相连接的一些值，弄成一个向量
  group_mutated <- as.numeric(d_f[1,mutated_indices])
  group_unmutated <- as.numeric(d_f[1,-mutated_indices])
  #进行检验
  wilcox_test_result <- wilcox.test(group_mutated, group_unmutated)
  p <- wilcox_test_result$p.value
  #结果写入result
  log2_fc <- median(as.numeric(group_mutated))-median(as.numeric(group_unmutated))
  result[i] <- paste(p,log2_fc,sep = ":")
}
  }
  return(paste(result,collapse = ";"))
}else{return(NA)
}}

#Func-2 在Func1的基础上只对cnnormal的样本进行计算(由于cnv文件只有chr1-22，这里也只能算chr1-22的基因)
#返回突变组相对非突变组的表达p值
p_expr_cnnormal <- function(genes,mutated_samples,nc_rownames,nc_hugosymbol,cnvfile){
  #如果genes不是NA，开弄
  if(all(!is.na(genes))){
    genes <- strsplit(genes,";")[[1]]
    mutated_samples <- strsplit(mutated_samples,",")[[1]]
    result <- rep(NA,length(genes))
    #所有的比较限于cnnormal的样本中
    #基因需要在cnv文件中有才行
    if(genes %in% rownames(cnvfile)){
    cn_samples <- names(cnvfile)[cnvfile[genes, ] == 0]
    cn_mutated_samples <- intersect(cn_samples,mutated_samples)
    cn_unmutated_samples <- setdiff(cn_samples,mutated_samples)
    #两个样本集不能为空
    if(!(is_empty(cn_mutated_samples) | is_empty(cn_unmutated_samples))){
    # 获取在表达矩阵中对应的列索引
    mutated_indices <- which(names(normalized_counts) %in% cn_mutated_samples)
    unmutated_indices <- which(names(normalized_counts) %in% cn_unmutated_samples)
    for(i in 1:length(genes)){
      #首先匹配到归一化counts矩阵该gene对应的行
      d_f <- data.frame()
      #如果该基因是ENSG开头,与在rownames里即可
      if(grepl("^ENSG",genes[i])){
        d_f <- normalized_counts[grepl(genes[i],nc_rownames),]
      }else{ #如果是Hugo_Symbol与Hugo_Symbol比对，需完全匹配
        d_f <- normalized_counts[nc_hugosymbol==genes[i],]
      }
      #如果能匹配上
      if(nrow(d_f)==1){
        # 将表达量分组
        group_mutated <- as.numeric(d_f[1,mutated_indices])
        group_unmutated <- as.numeric(d_f[1,unmutated_indices])
        #进行检验
        wilcox_test_result <- wilcox.test(group_mutated, group_unmutated)
        p <- wilcox_test_result$p.value
        #结果写入result
        log2_fc <- median(as.numeric(group_mutated))-median(as.numeric(group_unmutated))
        result[i] <- paste(p,log2_fc,sep = ":")
      }
    }
    return(paste(result,collapse = ";"))
    }else{return(NA)}
    }else{return(NA)}}else{return(NA)}}

#Func-3 输入突变的样本集合(Tumor_Sample_Barcode,逗号相连的一串值),返回Logrank检验的值(有正负)
#暂时弃用了 
if(0){
p_surv <- function(clinical_data,mutated_samples){
  d_f <- clinical_data
  mutated_samples <- strsplit(mutated_samples,",")[[1]]
  #分组
  d_f$Mut_type <- "Unmutated"
  d_f[mutated_samples,"Mut_type"] <- "Mutated"
  #做OS的logrank检验
  logrank_test <- survdiff(Surv(OS,OS_E) ~ Mut_type, data = d_f)
  p_value_os <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
  #做RFS的logrank检验
  logrank_test <- survdiff(Surv(RFS,RFS_E) ~ Mut_type, data = d_f)
  p_value_rfs <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
  return(paste(c(p_value_os,p_value_rfs),collapse = ";"))
}}

#Func-4 绘图函数1 输入df_sorted，返回柱状+热图组合图
#暂时弃用了 
if(0){
draw_picture <- function(df){
#按突变样本数从大到小重排数据框
rownames(df_sorted) <- paste(df_sorted$id,df_sorted$Target_Gene,sep = ";")
mat <- as.matrix(df_sorted[,c("expr_p","os_p","rfs_p")])

#设置热图颜色
col_fun = circlize::colorRamp2(c(0, 1), c("red", "white"))
values <- c("CDS" = "#3b6291",
            "UTR5" = "#943c39",
            "UTR3" = "#779043",
            "SS_CDS" = "#bf7334",
            "Promoter_CDS" = "#624c7c",
            "LncRNA" = "#388498",
            "SS_LncRNA" = "#82B0D2",
            "Promoter_LncRNA" = "#BEB8DC",
            "SmallRNA" = "#E7DAD2",
            "miRNA" = "#8ECFC9",
            "Enhancer" = "#FFBE7A",
            "CTCFbs" = "#FA7F6F",
            "Others" = "#B0B0B0")
df_sorted <- df_sorted %>%
  mutate(color = values[Type])
#设置中间列名
row_labels_mid <- sapply(rownames(df_sorted),function(x) strsplit(x,"[_;]")[[1]][3])
#设置行注释
row_ha <- rowAnnotation(
  bar = anno_barplot(df_sorted$mut.n,
                     gp = gpar(fill = df_sorted$color,
                               color="white"),
                     width=unit(10, "cm"),
                     axis_param=list(labels_rot=0)
  ),
  show_annotation_name = FALSE,
  Target_gene = anno_text(row_labels_mid,gp = gpar(fontsize = 6,fontface="bold")) # 大小
)
#设置最左边列名
row_labels <-  sapply(rownames(df_sorted),function(x) strsplit(x,"[_;]")[[1]][2])
#设置图例
#p_value
lgd1 = Legend(col_fun = col_fun,
              title = "p_value",
              labels_gp = gpar(col = "black"),
              title_gp = gpar(col = "black"))
#Type
lgd2 = Legend(direction = "vertical",
              title = "Type", 
              title_position = "topleft",
              labels = c("CDS",
                         "UTR5",
                         "UTR3",
                         "SS_CDS",
                         "Promoter_CDS",
                         "LncRNA",
                         "SS_LncRNA",
                         "Promoter_LncRNA",
                         "SmallRNA",
                         "miRNA",
                         "Enhancer",
                         "CTCFbs",
                         "Others"),
              grid_height = unit(0.4, "cm"), 
              grid_width = unit(4, "mm"),
              border = "black",
              legend_gp = gpar(fill = c("#3b6291",
                                                 "#943c39",
                                                 "#779043",
                                                 "#bf7334",
                                                 "#624c7c",
                                                 "#388498",
                                                 "#82B0D2",
                                                 "#BEB8DC",
                                                 "#E7DAD2",
                                                 "#8ECFC9",
                                                 "#FFBE7A",
                                                 "#FA7F6F",
                                                 "#B0B0B0")),
                                                 labels_gp = gpar(col = "black"),
              title_gp = gpar(col = "black")
)
pd = packLegend(lgd1,lgd2,
                direction = "vertical", # "vertical", "horizontal"
                #max_height = unit(20, "cm"), # 整个图例的最大高度
                #column_gap = unit(5, "mm"), #列间隔
                row_gap = unit(5, "mm")) #行间隔
heatmap <- Heatmap(mat,
                   col = col_fun,
                   rect_gp = gpar(col = "white", lwd = 3),
                   left_annotation = row_ha,
                   cluster_rows = F,
                   cluster_columns = F,
                   width = unit(2, "cm"),
                   row_names_side = "left", 
                   row_names_gp = gpar(fontsize = 6,fontface="bold"),
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 6,fontface="bold"),
                   column_names_rot = 0,#列名不旋转
                   column_names_centered = T,#列名居中
                   row_labels = row_labels,#设置行名
                   show_heatmap_legend = FALSE #不展示自动生成的热图图例
)
draw(heatmap,annotation_legend_list=pd)}}

#Func-5 绘图函数 输入df_sorted，返回图(q值)
draw_picture_q <- function(df){
  #按突变样本数从大到小重排数据框
  rownames(df_sorted) <- paste(df_sorted$id,df_sorted$Target_Gene,sep = ";")
  mat <- as.matrix(df_sorted[,c("cnnormal_expr_q")])
  
  #设置热图颜色
  col_fun = circlize::colorRamp2(c(0,1), c("red","white"))
  values <- c("CDS" = "#3b6291",
              "UTR5" = "#943c39",
              "UTR3" = "#779043",
              "SS_CDS" = "#bf7334",
              "Promoter_CDS" = "#624c7c",
              "LncRNA" = "#388498",
              "SS_LncRNA" = "#82B0D2",
              "Promoter_LncRNA" = "#BEB8DC",
              "SmallRNA" = "#E7DAD2",
              "miRNA" = "#8ECFC9",
              "Enhancer" = "#FFBE7A",
              "CTCFbs" = "#FA7F6F",
              "Others" = "#B0B0B0")
  df_sorted <- df_sorted %>%
    mutate(color = values[Type])
  #设置中间列名
  row_labels_mid <- sapply(rownames(df_sorted),function(x) strsplit(x,"[_;]")[[1]][3])
  #设置行注释
  row_ha <- rowAnnotation(
    bar = anno_barplot(df_sorted$mut.n,
                       gp = gpar(fill = df_sorted$color,
                                 color="white"),
                       width=unit(10, "cm"),
                       axis_param=list(labels_rot=0)
    ),
    show_annotation_name = FALSE,
    Target_gene = anno_text(row_labels_mid,gp = gpar(fontsize = 6,fontface="bold")) # 大小
  )
  #设置最左边列名
  row_labels <-  sapply(rownames(df_sorted),function(x) strsplit(x,"[_;]")[[1]][2])
  #设置图例
  #p_value
  lgd1 = Legend(col_fun = col_fun,
                title = "q_value",
                labels_gp = gpar(col = "black"),
                title_gp = gpar(col = "black"),
  )
  #Type
  lgd2 = Legend(direction = "vertical",
                title = "Type", 
                title_position = "topleft",
                labels = c("CDS",
                           "UTR5",
                           "UTR3",
                           "SS_CDS",
                           "Promoter_CDS",
                           "LncRNA",
                           "SS_LncRNA",
                           "Promoter_LncRNA",
                           "SmallRNA",
                           "miRNA",
                           "Enhancer",
                           "CTCFbs",
                           "Others"),
                grid_height = unit(0.4, "cm"), 
                grid_width = unit(4, "mm"),
                border = "black",
                legend_gp = gpar(fill = c("#3b6291",
                                                   "#943c39",
                                                   "#779043",
                                                   "#bf7334",
                                                   "#624c7c",
                                                   "#388498",
                                                   "#82B0D2",
                                                   "#BEB8DC",
                                                   "#E7DAD2",
                                                   "#8ECFC9",
                                                   "#FFBE7A",
                                                   "#FA7F6F",
                                                   "#B0B0B0")),
                                                   labels_gp = gpar(col = "black"),
                title_gp = gpar(col = "black")
  )
  pd = packLegend(lgd1,lgd2,
                  direction = "vertical", # "vertical", "horizontal"
                  #max_height = unit(20, "cm"), # 整个图例的最大高度
                  #column_gap = unit(5, "mm"), #列间隔
                  row_gap = unit(5, "mm")) #行间隔
  heatmap <- Heatmap(mat,
                     na_col = "grey",
                     col = col_fun,
                     rect_gp = gpar(col = "white", lwd = 3),
                     left_annotation = row_ha,
                     cluster_rows = F,
                     cluster_columns = F,
                     width = unit(1, "cm"),
                     row_names_side = "left", 
                     row_names_gp = gpar(fontsize = 6,fontface="bold"),
                     column_names_side = "top",
                     column_names_gp = gpar(fontsize = 10,fontface="bold"),
                     column_names_rot = 0,#列名不旋转
                     column_names_centered = T,#列名居中
                     row_labels = row_labels,#设置行名
                     column_labels = c("expr"),#设置列名
                     show_heatmap_legend = FALSE #不展示自动生成的热图图例
  )
  draw(heatmap,annotation_legend_list=pd)}
}

#===============================================================================

#1 导入和预处理文件
if(1){
#得到基因坐标文件
all_gene_bed_hg38 <- read.delim(all_gene_bed_hg38_file_path,sep=" ",header=F)
names(all_gene_bed_hg38) <- c("chr","start","end","gene")
#得到config文件
config <- read.delim(config_filepath,header = F)
config$V1 <- sapply(config$V1, function(x) {str_remove(x, "(LC).*")})
#得到clinical_data文件
clinical_data <- read.delim(clinical_data_file_path)
sample_id_986 <- paste(config$V1,"LC",sep = "")
clinical_data <- clinical_data[clinical_data$Sample_ID %in% sample_id_986,]
rownames(clinical_data) <- clinical_data$Sample_ID
#得到肺癌驱动突变
lung_drivers_df <- read.csv(lung_drivers_df_filepath)
lung_drivers <- lung_drivers_df$Gene_Symbol

#得到所有癌的驱动突变
oncokb_cancergene_list_df <- read.delim(oncokb_cancergene_list_filepath)
oncokb_cancergene_list <- oncokb_cancergene_list_df[,1]
cancer_gene_737_df <- read.csv(cancer_gene_737_filepath)
cancer_gene_737 <- cancer_gene_737_df[,1]
cosmic_census_cancer_gene_df <- read.delim(cosmic_census_cancer_gene_filepath)
cosmic_census_cancer_gene <- cosmic_census_cancer_gene_df[,c("Gene.Symbol")]
all_cancer_gene <- union(cancer_gene_737,cosmic_census_cancer_gene)
all_cancer_gene <- union(all_cancer_gene,oncokb_cancergene_list)
all_cancer_gene <- union(all_cancer_gene,lung_drivers)
#导入注释文件 by THT
annotations_hg38 <- fread(annotations_hg38_filepath)
#设置成Granges对象(因为bed是左闭右开，所以结尾减去1)
gr_anno_hg38 <- GRanges(seqnames = annotations_hg38$V1,ranges = IRanges(start = annotations_hg38$V2,end = annotations_hg38$V3-1))

#得到TAD文件及其GRanges对象
TAD_withgenes_hg38 <- fread(TAD_withgenes_hg38_filepath)
gr_TAD_withgenes_hg38 <- GRanges(seqnames = TAD_withgenes_hg38$chr,ranges = IRanges(start=TAD_withgenes_hg38$start,end=TAD_withgenes_hg38$end -1))

#得到归一化后的counts
normalized_counts <- fread(normalized_counts_filepath,data.table = F)
#fread不支持行名处理一下
rownames(normalized_counts) <- normalized_counts$V1
normalized_counts$V1 <- NULL

#只留下config有的样本Tumor的列
mid <- paste0(config$V1,"_T")
pattern <- paste(mid,collapse = "|")
normalized_counts <- normalized_counts[,grepl(pattern,colnames(normalized_counts))]
#保留三位小数
#normalized_counts <- as.data.frame(lapply(normalized_counts, round, 3))
#改一下列名
colnames(normalized_counts) <- sapply(colnames(normalized_counts),function(x) paste0(strsplit(x,"_")[[1]][3],strsplit(x,"_")[[1]][2]))

#RAS/RAF/MAPK通路 相关基因
RPA_genes <- c("KRAS","EGFR","BRAF","ERBB2", 
               "MET", "RIT1", "NRAS", "RAF1", 
               "HRAS", "ARAF", "MAP2K1", "SOS1",
               "NF1","RASA1","ALK","ROS1","RET", 
               "MET","NTRK2","FGFR1" ,"MAPK1")
#拷贝数结果
cnv <- fread(cnv_filepath,data.table = F) #里面的值是log2(CN)-1
cnv_1 <- cnv
rownames(cnv_1) <- cnv_1$`Gene Symbol`
cnv_1[,1:3] <- NULL

#每个点突变的突变特征
all_singlemuts_sigs <- fread(all_singlemuts_sigs_filepath,data.table = F)
all_singlemuts_sigs <- all_singlemuts_sigs %>%
  select(Sample.Names,Chr,Pos,MutationType,Major_SBS) %>%
  mutate(Chr = paste0("chr",Chr))

#Mappability及IR
GEM_mappability <- fread(GEM_mappability_filepath) #0-based 左开右闭
names(GEM_mappability) <- c("chr","start","end")
encode_blacklist_regions <- fread(Encode_blacklist_regions_filepath) #0-based 左开右闭
names(encode_blacklist_regions) <- c("chr","start","end")
IR <- readRDS(IR_filepath) #1-based 双闭
IR <- reduce(IR)
}

#===============================================================================

#2 得到输出数据框
#2-1 处理enpmi数据框，得到表格1 enpm_tad_output
#表格1说明
#Location 例：chr1:110002
#Overlap_element Overlap的注释元件,除enhancer外一个点仅一个
#Type 包括Enhancer,"CDS","Promoter_CDS等
#mut.n 共有多少个样本有该突变
#Mutated_Samples 具体是哪些样本有该突变
#Mappable 布尔值，该突变是否位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#In_IR 布尔值，该突变是否位于一段回文序列里
#Target_Gene 注释元件的靶基因,others对应的靶基因是相同tad内基因
#raw_expr_p 突变组与非突变组表达值wilcoxin检验的p值，无拷贝数矫正
#raw_expr_log2fc 突变组与非突变组表达值的log2FC，无拷贝数矫正
#cnnormal_expr_p 突变组与非突变组表达值wilcoxin检验的p值，有拷贝数矫正
#cnnormal_expr_log2fc 突变组与非突变组表达值的log2FC，有拷贝数矫正
#raw_expr_q 突变组与非突变组表达值fdr q值，无拷贝数矫正
#cnnormal_expr_q 突变组与非突变组表达值fdr q值，有拷贝数矫正
if(1){
#2-1-1 导入vep注释文件，包含location
  vep_result <- fread(enpmi_filepath)
  vep_result$location <- paste(vep_result$V2,vep_result$V3,sep = ":")
  vep_result <- vep_result[,-c(2,3)]
  colnames(vep_result) <- c("Tumor_Sample_Barcode","ref","alt","Location")
#去掉indel
vep_result <- vep_result %>%
    filter(nchar(ref) == 1 & nchar(alt) == 1)

#2-1-2 添加Overlap_element
if(1){
  #作为query的注释文件已经弄成在前面弄成granges格式了
  #利用vep_result弄一个subject GRanges文件
  if(1){
  seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
  #vep来自vcf，1-based，要减去1
  start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))-1
  end <- start
  gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
  #找到hits
  hits <- findOverlaps(gr_anno_hg38,gr_vep_result)
  anno_rownumbers <- queryHits(hits)
  vep_rownumbers <- subjectHits(hits)
  #设置一些中间向量减少内存使用
  result <- rep(NA,2447)
  anno_elements <- annotations_hg38$V4
  for(i in 1:length(hits)){
    #该hit result的值
    a <- result[vep_rownumbers[i]]
    #该hit 对应的anno的第四列的值
    b <- anno_elements[anno_rownumbers[i]]
    #如果还没赋过值，直接赋值
    if(is.na(a)){
      result[vep_rownumbers[i]] <- b}
    #否则加在后面
    else{
      result[vep_rownumbers[i]] <- paste(a,b,sep = " ; ")}
  }
  vep_result$Overlap_element <- result
#Overlap_element可以有多个，以;相分隔，将其拆开成每个注释一行
vep_result <- vep_result %>%
  separate_rows(Overlap_element, sep = ";")
}
#2-1-3 添加一列Type
if(1){
vep_result <- vep_result %>%
  mutate(Type = case_when(
    str_detect(Overlap_element, "gc46_pc.cds") ~ "CDS",
    str_detect(Overlap_element, "gc46_pc.ss_cds") ~ "SS_CDS",
    str_detect(Overlap_element, "gc46_pc.5utr") ~ "UTR5",
    str_detect(Overlap_element, "gc46_pc.3utr") ~ "UTR3",
    str_detect(Overlap_element, "gc46.promoter") ~ "Promoter_CDS",
    str_detect(Overlap_element, "lncrna_exons") ~ "LncRNA",
    str_detect(Overlap_element, "ss_lncrna") ~ "SS_LncRNA",
    str_detect(Overlap_element, "lncrna.promoter") ~ "Promoter_LncRNA",
    str_detect(Overlap_element, "smallrna") ~ "SmallRNA",
    str_detect(Overlap_element, "mirna") ~ "miRNA",
    str_detect(Overlap_element, "enhancers") ~ "Enhancer",
    str_detect(Overlap_element, "ctcfbs") ~ "CTCFbs",
    TRUE ~ "Others"  # 其他情况赋值为 Others
  ))  
}

#2-1-4 取出有用的几列
vep_result <- vep_result[,c("Tumor_Sample_Barcode","Location","Overlap_element","Type")]

#2-1-5 把相同位置且Overlap_element相同的突变整合成一行，添加新列 mut.n(有多少个样本含该突变)和Mutated_Samples(具体是哪些样本)
vep_result <- vep_result %>%
  group_by(across(-Tumor_Sample_Barcode)) %>%
  mutate(mut.n=n(),Mutated_Samples=paste(Tumor_Sample_Barcode,collapse = ","))  %>%
  group_by(across(-Tumor_Sample_Barcode)) %>%
  summarise()

#2-1-6 添加Target_gene(Hugo_Symbol或ENSG)
#首先添加除Enhancer外的靶基因
matching_types <- c("CDS","UTR3","LncRNA", "miRNA", "Promoter_CDS", "Promoter_LncRNA", 
                    "SmallRNA", "SS_CDS", "SS_LncRNA", "UTR5")
vep_result$Target_Gene <- ifelse(vep_result$Type %in% matching_types,sapply(strsplit(vep_result$Overlap_element, "::"), `[`, 3),NA)
#然后添加Enhancer的
mid <- vep_result[vep_result$Type == "Enhancer",]
vep_result$Target_Gene[vep_result$Type == "Enhancer"] <-sapply(strsplit(mid$Overlap_element, "::"), `[`, 4)

#最后添加Others的，Target_gene取同一个TAD的基因
if(1){
df1 <- vep_result%>%
  mutate(chr = sapply(Location, function(x) strsplit(x, "[:]")[[1]][1]),
         start = sapply(Location, function(x) as.integer(strsplit(x, "[:]")[[1]][2])-1),
         end = sapply(Location, function(x) as.integer(strsplit(x, "[:]")[[1]][2])-1),
  )
#把Others对应的行取出来
Others <- df1[df1$Type == "Others",]
#做Overlap
gr_others <- GRanges(seqnames = Others$chr,ranges = IRanges(start =Others$start ,end = Others$end))
hits <- findOverlaps(gr_others,gr_TAD_withgenes_hg38)
TAD_rownumber <- subjectHits(hits)
other_rownumber <- queryHits(hits)
#得到Target_gene
genes <- TAD_withgenes_hg38$genes_within
Others$Target_Gene[other_rownumber] <- genes[TAD_rownumber]
vep_result$Target_Gene[vep_result$Type == "Others"] <- Others$Target_Gene}
#Others可以有多个Target_gene,拆开成一个一行
vep_result <- vep_result %>%
  separate_rows(Target_Gene, sep = ";")
#去掉mut.n小于10的行
vep_result <- vep_result[vep_result$mut.n >=10,]
#有些tad里有重复的基因，去除
vep_result <- vep_result[!duplicated(vep_result),]

#2-1-7添加 Mappable
#2-1-7-1得到gr_mappable(1-based)
if(1){
#是0-based 左闭右开 start+1 end+1-1
gr_gem_mappabilty <- GRanges(seqnames = GEM_mappability$chr,ranges = IRanges(start = GEM_mappability$start+1,end = GEM_mappability$end))
gr_gem_mappabilty <- GenomicRanges::reduce(gr_gem_mappabilty)
#是0-based 左闭右开 start+1 end+1-1
gr_encode_blacklist <- GRanges(seqnames = encode_blacklist_regions$chr,ranges = IRanges(start = encode_blacklist_regions$start+1,end = encode_blacklist_regions$end))
gr_encode_blacklist <- GenomicRanges::reduce(gr_encode_blacklist)
#用gem减去encode_blacklist
gr_mappable <- GenomicRanges::setdiff(gr_gem_mappabilty, gr_encode_blacklist)}
#2-1-7-2得到gr_vep_result(和2-2中的不一样，这里是1-based)
if(1){
  seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
  #vep来自vcf，1-based
  start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
  end <- start
  gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
hits_map_vepresult <- findOverlaps(gr_mappable,gr_vep_result)
vep_result[,"Mappable"] <- FALSE
vep_result[subjectHits(hits_map_vepresult),"Mappable"] <- TRUE

#2-1-8添加 In_IR
#2-1-8-1 IR本身就是GRanges对象了，不用获得
#2-1-8-2得到gr_vep_result(和2-2中的不一样，这里是1-based)
if(1){
  seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
  #vep来自vcf，1-based
  start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
  end <- start
  gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
hits_IR_vep <- findOverlaps(IR,gr_vep_result)
vep_result[,"IR"] <- FALSE
vep_result[subjectHits(hits_IR_vep),"IR"] <- TRUE

#2-1-9 添加靶基因expr的p和log2fc值(without cnv correction) #1h
nc_rownames <- rownames(normalized_counts)
nc_hugosymbol <- sapply(nc_rownames,function(x) strsplit(x,"\\|")[[1]][2])
#result <- sapply(1:nrow(vep_result),function(i) p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol))
cl <- makeCluster(20)
registerDoParallel(cl)
result <- foreach(i = 1:nrow(vep_result), .combine = c) %dopar% {
  p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol)
}
stopCluster(cl)
vep_result$expr_p_log2fc <- result
vep_result <- vep_result %>%
  separate(expr_p_log2fc,into = c("raw_expr_p","raw_expr_log2fc"),sep=":")
  
#2-1-10 添加靶基因expr的p_cnnormal和log2fc值(with cnv correction) #1h
nc_rownames <- rownames(normalized_counts)
nc_hugosymbol <- sapply(nc_rownames,function(x) strsplit(x,"\\|")[[1]][2])
#result <- sapply(1:nrow(vep_result),function(i) p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol))
cl <- makeCluster(20)
registerDoParallel(cl)
result <- foreach(i = 1:nrow(vep_result), .combine = c) %dopar% {
  library(rlang)
  p_expr_cnnormal(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol,cnv_1)
}
stopCluster(cl)
vep_result$expr_p_log2fc_cnnormal <- result
vep_result1 <- vep_result %>%
  separate(expr_p_log2fc_cnnormal,into = c("cnnormal_expr_p","cnnormal_expr_log2fc"),sep=":")
#2-1-11 由于有一些基因在超过50%的基因中表达都为0，需要把这些基因算出来的expr_p去掉
#得到这些基因在表达矩阵中对应的行名
# 计算每行值为0的列的数量
zero_count <- apply(normalized_counts == 0, 1, sum)
# 计算总列数
total_columns <- ncol(normalized_counts)
# 筛选出小于50%列值为0的行的行名
rownames  <- rownames(normalized_counts[zero_count < total_columns / 2, ])
# 使用 strsplit 按照 "|" 分割
split_vec <- strsplit(rownames, "\\|")
# 扁平化为一个向量
flattened_vec <- unlist(split_vec)
#去掉vep_result1中的这些基因的expr
vep_result1 <- vep_result1 %>%
  mutate(raw_expr_p=ifelse(Target_Gene %in% flattened_vec,raw_expr_p,NA)) %>%
  mutate(raw_expr_log2fc=ifelse(Target_Gene %in% flattened_vec,raw_expr_log2fc,NA)) %>%
  mutate(cnnormal_expr_p=ifelse(Target_Gene %in% flattened_vec,cnnormal_expr_p,NA)) %>%
  mutate(cnnormal_expr_log2fc=ifelse(Target_Gene %in% flattened_vec,cnnormal_expr_log2fc,NA))
#2-1-12 计算raw expr q值
#raw_expr_p的NA可能是因为Func1有点问题，实际上是字符串
vep_result1$raw_expr_p <- ifelse(!(vep_result1$raw_expr_p=="NA"),vep_result1$raw_expr_p,NA)
vep_result1$raw_expr_p <- as.numeric(vep_result1$raw_expr_p)
vep_result1$raw_expr_q <- p.adjust(vep_result1$raw_expr_p,method = "BH")
#2-1-13 计算cnnormal expr q值
vep_result1$cnnormal_expr_p <- as.numeric(vep_result1$cnnormal_expr_p)
vep_result1$cnnormal_expr_q <- p.adjust(vep_result1$cnnormal_expr_p,method = "BH")
enpm_output_tad <- vep_result1
write.table(enpm_output_tad,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_output_tad.tsv",row.names = F,quote = F,sep = "\t")
}

#2-2 处理enpmi数据框，得到表格2 enpm_5M_output
#表格2说明
#Location 例：chr1:110002
#Overlap_element Overlap的注释元件,除enhancer外一个点仅一个
#Type 包括Enhancer,"CDS","Promoter_CDS等
#mut.n 共有多少个样本有该突变
#Mutated_Samples 具体是哪些样本有该突变
#Mappable 布尔值，该突变是否位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#In_IR 布尔值，该突变是否位于一段回文序列里
#Target_Gene 注释元件的靶基因,others对应的靶基因是距离5M内基因
#raw_expr_p 突变组与非突变组表达值wilcoxin检验的p值，无拷贝数矫正
#raw_expr_log2fc 突变组与非突变组表达值的log2FC，无拷贝数矫正
#cnnormal_expr_p 突变组与非突变组表达值wilcoxin检验的p值，有拷贝数矫正
#cnnormal_expr_log2fc 突变组与非突变组表达值的log2FC，有拷贝数矫正
#raw_expr_q 突变组与非突变组表达值fdr q值，无拷贝数矫正
#cnnormal_expr_q 突变组与非突变组表达值fdr q值，有拷贝数矫正
if(1){
  #2-2-1 导入vep注释文件，包含location
  vep_result <- fread(enpmi_filepath)
  vep_result$location <- paste(vep_result$V2,vep_result$V3,sep = ":")
  vep_result <- vep_result[,-c(2,3)]
  colnames(vep_result) <- c("Tumor_Sample_Barcode","ref","alt","Location")
  #去掉indel
  vep_result <- vep_result %>%
    filter(nchar(ref) == 1 & nchar(alt) == 1)
  
  #2-2-2 添加Overlap_element
  if(1){
    #作为query的注释文件已经弄成在前面弄成granges格式了
    #利用vep_result弄一个subject GRanges文件
    if(1){
      seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
      #vep来自vcf，1-based，要减去1
      start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))-1
      end <- start
      gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
    #找到hits
    hits <- findOverlaps(gr_anno_hg38,gr_vep_result)
    anno_rownumbers <- queryHits(hits)
    vep_rownumbers <- subjectHits(hits)
    #设置一些中间向量减少内存使用
    result <- rep(NA,2447)
    anno_elements <- annotations_hg38$V4
    for(i in 1:length(hits)){
      #该hit result的值
      a <- result[vep_rownumbers[i]]
      #该hit 对应的anno的第四列的值
      b <- anno_elements[anno_rownumbers[i]]
      #如果还没赋过值，直接赋值
      if(is.na(a)){
        result[vep_rownumbers[i]] <- b}
      #否则加在后面
      else{
        result[vep_rownumbers[i]] <- paste(a,b,sep = " ; ")}
    }
    vep_result$Overlap_element <- result
    #Overlap_element可以有多个，以;相分隔，将其拆开成每个注释一行
    vep_result <- vep_result %>%
      separate_rows(Overlap_element, sep = ";")
  }
  #2-2-3 添加一列Type
  if(1){
    vep_result <- vep_result %>%
      mutate(Type = case_when(
        str_detect(Overlap_element, "gc46_pc.cds") ~ "CDS",
        str_detect(Overlap_element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(Overlap_element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(Overlap_element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(Overlap_element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(Overlap_element, "lncrna_exons") ~ "LncRNA",
        str_detect(Overlap_element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(Overlap_element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(Overlap_element, "smallrna") ~ "SmallRNA",
        str_detect(Overlap_element, "mirna") ~ "miRNA",
        str_detect(Overlap_element, "enhancers") ~ "Enhancer",
        str_detect(Overlap_element, "ctcfbs") ~ "CTCFbs",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      ))  
  }
  
  #2-2-4 取出有用的几列
  vep_result <- vep_result[,c("Tumor_Sample_Barcode","Location","Overlap_element","Type")]
  
  #2-2-5 把相同位置且Overlap_element相同的突变整合成一行，添加新列 mut.n(有多少个样本含该突变)和Mutated_Samples(具体是哪些样本)
  vep_result <- vep_result %>%
    group_by(across(-Tumor_Sample_Barcode)) %>%
    mutate(mut.n=n(),Mutated_Samples=paste(Tumor_Sample_Barcode,collapse = ","))  %>%
    group_by(across(-Tumor_Sample_Barcode)) %>%
    summarise()
  
  #2-2-6 添加Target_gene(Hugo_Symbol或ENSG)
  #首先添加除Enhancer外的靶基因
  matching_types <- c("CDS","UTR3","LncRNA", "miRNA", "Promoter_CDS", "Promoter_LncRNA", 
                      "SmallRNA", "SS_CDS", "SS_LncRNA", "UTR5")
  vep_result$Target_Gene <- ifelse(vep_result$Type %in% matching_types,sapply(strsplit(vep_result$Overlap_element, "::"), `[`, 3),NA)
  #然后添加Enhancer的
  mid <- vep_result[vep_result$Type == "Enhancer",]
  vep_result$Target_Gene[vep_result$Type == "Enhancer"] <-sapply(strsplit(mid$Overlap_element, "::"), `[`, 4)
  
  #最后添加Others的，Target_gene取距离5M内的基因
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_gene <- GRanges(seqnames = all_gene_bed_hg38$chr,ranges = IRanges(start = all_gene_bed_hg38$start+1,end = all_gene_bed_hg38$end))
    #1-based 左闭右闭，左不变，右不变
    if(1){
      seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
      #vep来自vcf，1-based，要减去1
      start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))-1
      end <- start
      gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
    hits_gene_vep <- findOverlaps(query = gr_gene,subject = gr_vep_result,maxgap = 5000000)
    #添加Others的target_gene
    #每个Others对应多个target_gene,先把它们整理成每个id一行
    vep_result$id <- paste(rownames(vep_result), vep_result$Location, sep = "_")
    others_target <- data.frame(id=vep_result$id[subjectHits(hits_gene_vep)],others_target=all_gene_bed_hg38$gene[queryHits(hits_gene_vep)])
    others_target <- others_target %>%
      group_by(id) %>%
      summarise(others_targets=paste(others_target,collapse = ";"))
    vep_result <- vep_result %>%
      left_join(others_target,by="id") %>%
      mutate(Target_Gene=ifelse(is.na(Target_Gene),others_targets,Target_Gene)) %>%
      separate_rows(Target_Gene,sep = ";") %>%
      select(-id,-others_targets)
  }
  #去掉mut.n小于10的行
  vep_result <- vep_result[vep_result$mut.n >=10,]
  #有些tad里有重复的基因，去除
  vep_result <- vep_result[!duplicated(vep_result),]
  
  #2-2-7添加 Mappable
  #2-2-7-1得到gr_mappable(1-based)
  if(1){
    #是0-based 左闭右开 start+1 end+1-1
    gr_gem_mappabilty <- GRanges(seqnames = GEM_mappability$chr,ranges = IRanges(start = GEM_mappability$start+1,end = GEM_mappability$end))
    gr_gem_mappabilty <- GenomicRanges::reduce(gr_gem_mappabilty)
    #是0-based 左闭右开 start+1 end+1-1
    gr_encode_blacklist <- GRanges(seqnames = encode_blacklist_regions$chr,ranges = IRanges(start = encode_blacklist_regions$start+1,end = encode_blacklist_regions$end))
    gr_encode_blacklist <- GenomicRanges::reduce(gr_encode_blacklist)
    #用gem减去encode_blacklist
    gr_mappable <- GenomicRanges::setdiff(gr_gem_mappabilty, gr_encode_blacklist)}
  #2-1-7-2得到gr_vep_result(和2-2中的不一样，这里是1-based)
  if(1){
    seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
    #vep来自vcf，1-based
    start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
    end <- start
    gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
  hits_map_vepresult <- findOverlaps(gr_mappable,gr_vep_result)
  vep_result[,"Mappable"] <- FALSE
  vep_result[subjectHits(hits_map_vepresult),"Mappable"] <- TRUE
  
  #2-2-8添加 In_IR
  #2-2-8-1 IR本身就是GRanges对象了，不用获得
  #2-2-8-2得到gr_vep_result(和2-2中的不一样，这里是1-based)
  if(1){
    seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
    #vep来自vcf，1-based
    start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
    end <- start
    gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
  hits_IR_vep <- findOverlaps(IR,gr_vep_result)
  vep_result[,"IR"] <- FALSE
  vep_result[subjectHits(hits_IR_vep),"IR"] <- TRUE
  
  #2-2-9 添加靶基因expr的p和log2fc值(without cnv correction) #1h
  nc_rownames <- rownames(normalized_counts)
  nc_hugosymbol <- sapply(nc_rownames,function(x) strsplit(x,"\\|")[[1]][2])
  #result <- sapply(1:nrow(vep_result),function(i) p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol))
  cl <- makeCluster(40)
  registerDoParallel(cl)
  result <- foreach(i = 1:nrow(vep_result), .combine = c) %dopar% {
    p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol)
  }
  stopCluster(cl)
  vep_result$expr_p_log2fc <- result
  vep_result <- vep_result %>%
    separate(expr_p_log2fc,into = c("raw_expr_p","raw_expr_log2fc"),sep=":")
    
#2-2-10 添加靶基因expr的p_cnnormal和log2fc值(with cnv correction) #1h
    nc_rownames <- rownames(normalized_counts)
  nc_hugosymbol <- sapply(nc_rownames,function(x) strsplit(x,"\\|")[[1]][2])
  #result <- sapply(1:nrow(vep_result),function(i) p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol))
  cl <- makeCluster(40)
  registerDoParallel(cl)
  result <- foreach(i = 1:nrow(vep_result), .combine = c) %dopar% {
    library(rlang)
    p_expr_cnnormal(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol,cnv_1)
  }
  stopCluster(cl)
  vep_result$expr_p_log2fc_cnnormal <- result
  vep_result1 <- vep_result %>%
    separate(expr_p_log2fc_cnnormal,into = c("cnnormal_expr_p","cnnormal_expr_log2fc"),sep=":")
#2-2-11 由于有一些基因在超过50%的基因中表达都为0，需要把这些基因算出来的expr_p去掉
  #得到这些基因在表达矩阵中对应的行名
  # 计算每行值为0的列的数量
  zero_count <- apply(normalized_counts == 0, 1, sum)
  # 计算总列数
  total_columns <- ncol(normalized_counts)
  # 筛选出小于50%列值为0的行的行名
  rownames  <- rownames(normalized_counts[zero_count < total_columns / 2, ])
  # 使用 strsplit 按照 "|" 分割
  split_vec <- strsplit(rownames, "\\|")
  # 扁平化为一个向量
  flattened_vec <- unlist(split_vec)
  #去掉vep_result1中的这些基因的expr
  vep_result1 <- vep_result1 %>%
    mutate(raw_expr_p=ifelse(Target_Gene %in% flattened_vec,raw_expr_p,NA)) %>%
    mutate(raw_expr_log2fc=ifelse(Target_Gene %in% flattened_vec,raw_expr_log2fc,NA)) %>%
    mutate(cnnormal_expr_p=ifelse(Target_Gene %in% flattened_vec,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_log2fc=ifelse(Target_Gene %in% flattened_vec,cnnormal_expr_log2fc,NA))
#2-2-12 计算raw expr q值
  #raw_expr_p的NA可能是因为Func1有点问题，实际上是字符串
  vep_result1$raw_expr_p <- ifelse(!(vep_result1$raw_expr_p=="NA"),vep_result1$raw_expr_p,NA)
  vep_result1$raw_expr_p <- as.numeric(vep_result1$raw_expr_p)
  vep_result1$raw_expr_q <- p.adjust(vep_result1$raw_expr_p,method = "BH")
#2-2-13 计算cnnormal expr q值
  vep_result1$cnnormal_expr_p <- as.numeric(vep_result1$cnnormal_expr_p)
  vep_result1$cnnormal_expr_q <- p.adjust(vep_result1$cnnormal_expr_p,method = "BH")
  enpm_output_5M <- vep_result1
  write.table(enpm_output_5M,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_output_5M.tsv",row.names = F,quote = F,sep = "\t")
}

#2-3 处理enpmi数据框，得到表格2 enpm_10M_output
#表格2说明
#Location 例：chr1:110002
#Overlap_element Overlap的注释元件,除enhancer外一个点仅一个
#Type 包括Enhancer,"CDS","Promoter_CDS等
#mut.n 共有多少个样本有该突变
#Mutated_Samples 具体是哪些样本有该突变
#Mappable 布尔值，该突变是否位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#In_IR 布尔值，该突变是否位于一段回文序列里
#Target_Gene 注释元件的靶基因,others对应的靶基因是距离10M内基因
#raw_expr_p 突变组与非突变组表达值wilcoxin检验的p值，无拷贝数矫正
#raw_expr_log2fc 突变组与非突变组表达值的log2FC，无拷贝数矫正
#cnnormal_expr_p 突变组与非突变组表达值wilcoxin检验的p值，有拷贝数矫正
#cnnormal_expr_log2fc 突变组与非突变组表达值的log2FC，有拷贝数矫正
#raw_expr_q 突变组与非突变组表达值fdr q值，无拷贝数矫正
#cnnormal_expr_q 突变组与非突变组表达值fdr q值，有拷贝数矫正
if(1){
  #2-3-1 导入vep注释文件，包含location
  vep_result <- fread(enpmi_filepath)
  vep_result$location <- paste(vep_result$V2,vep_result$V3,sep = ":")
  vep_result <- vep_result[,-c(2,3)]
  colnames(vep_result) <- c("Tumor_Sample_Barcode","ref","alt","Location")
  #去掉indel
  vep_result <- vep_result %>%
    filter(nchar(ref) == 1 & nchar(alt) == 1)
  
  #2-3-2 添加Overlap_element
  if(1){
    #作为query的注释文件已经弄成在前面弄成granges格式了
    #利用vep_result弄一个subject GRanges文件
    if(1){
      seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
      #vep来自vcf，1-based，要减去1
      start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))-1
      end <- start
      gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
    #找到hits
    hits <- findOverlaps(gr_anno_hg38,gr_vep_result)
    anno_rownumbers <- queryHits(hits)
    vep_rownumbers <- subjectHits(hits)
    #设置一些中间向量减少内存使用
    result <- rep(NA,2447)
    anno_elements <- annotations_hg38$V4
    for(i in 1:length(hits)){
      #该hit result的值
      a <- result[vep_rownumbers[i]]
      #该hit 对应的anno的第四列的值
      b <- anno_elements[anno_rownumbers[i]]
      #如果还没赋过值，直接赋值
      if(is.na(a)){
        result[vep_rownumbers[i]] <- b}
      #否则加在后面
      else{
        result[vep_rownumbers[i]] <- paste(a,b,sep = " ; ")}
    }
    vep_result$Overlap_element <- result
    #Overlap_element可以有多个，以;相分隔，将其拆开成每个注释一行
    vep_result <- vep_result %>%
      separate_rows(Overlap_element, sep = ";")
  }
  #2-3-3 添加一列Type
  if(1){
    vep_result <- vep_result %>%
      mutate(Type = case_when(
        str_detect(Overlap_element, "gc46_pc.cds") ~ "CDS",
        str_detect(Overlap_element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(Overlap_element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(Overlap_element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(Overlap_element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(Overlap_element, "lncrna_exons") ~ "LncRNA",
        str_detect(Overlap_element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(Overlap_element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(Overlap_element, "smallrna") ~ "SmallRNA",
        str_detect(Overlap_element, "mirna") ~ "miRNA",
        str_detect(Overlap_element, "enhancers") ~ "Enhancer",
        str_detect(Overlap_element, "ctcfbs") ~ "CTCFbs",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      ))  
  }
  
  #2-3-4 取出有用的几列
  vep_result <- vep_result[,c("Tumor_Sample_Barcode","Location","Overlap_element","Type")]
  
  #2-3-5 把相同位置且Overlap_element相同的突变整合成一行，添加新列 mut.n(有多少个样本含该突变)和Mutated_Samples(具体是哪些样本)
  vep_result <- vep_result %>%
    group_by(across(-Tumor_Sample_Barcode)) %>%
    mutate(mut.n=n(),Mutated_Samples=paste(Tumor_Sample_Barcode,collapse = ","))  %>%
    group_by(across(-Tumor_Sample_Barcode)) %>%
    summarise()
  
  #2-3-6 添加Target_gene(Hugo_Symbol或ENSG)
  #首先添加除Enhancer外的靶基因
  matching_types <- c("CDS","UTR3","LncRNA", "miRNA", "Promoter_CDS", "Promoter_LncRNA", 
                      "SmallRNA", "SS_CDS", "SS_LncRNA", "UTR5")
  vep_result$Target_Gene <- ifelse(vep_result$Type %in% matching_types,sapply(strsplit(vep_result$Overlap_element, "::"), `[`, 3),NA)
  #然后添加Enhancer的
  mid <- vep_result[vep_result$Type == "Enhancer",]
  vep_result$Target_Gene[vep_result$Type == "Enhancer"] <-sapply(strsplit(mid$Overlap_element, "::"), `[`, 4)
  
  #最后添加Others的，Target_gene取距离5M内的基因
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_gene <- GRanges(seqnames = all_gene_bed_hg38$chr,ranges = IRanges(start = all_gene_bed_hg38$start+1,end = all_gene_bed_hg38$end))
    #1-based 左闭右闭，左不变，右不变
    if(1){
      seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
      #vep来自vcf，1-based，要减去1
      start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))-1
      end <- start
      gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
    hits_gene_vep <- findOverlaps(query = gr_gene,subject = gr_vep_result,maxgap = 10000000)
    #添加Others的target_gene
    #每个Others对应多个target_gene,先把它们整理成每个id一行
    vep_result$id <- paste(rownames(vep_result), vep_result$Location, sep = "_")
    others_target <- data.frame(id=vep_result$id[subjectHits(hits_gene_vep)],others_target=all_gene_bed_hg38$gene[queryHits(hits_gene_vep)])
    others_target <- others_target %>%
      group_by(id) %>%
      summarise(others_targets=paste(others_target,collapse = ";"))
    vep_result <- vep_result %>%
      left_join(others_target,by="id") %>%
      mutate(Target_Gene=ifelse(is.na(Target_Gene),others_targets,Target_Gene)) %>%
      separate_rows(Target_Gene,sep = ";") %>%
      select(-id,-others_targets)
  }
  #去掉mut.n小于10的行
  vep_result <- vep_result[vep_result$mut.n >=10,]
  #有些tad里有重复的基因，去除
  vep_result <- vep_result[!duplicated(vep_result),]
  
  #2-3-7添加 Mappable
  #2-3-7-1得到gr_mappable(1-based)
  if(1){
    #是0-based 左闭右开 start+1 end+1-1
    gr_gem_mappabilty <- GRanges(seqnames = GEM_mappability$chr,ranges = IRanges(start = GEM_mappability$start+1,end = GEM_mappability$end))
    gr_gem_mappabilty <- GenomicRanges::reduce(gr_gem_mappabilty)
    #是0-based 左闭右开 start+1 end+1-1
    gr_encode_blacklist <- GRanges(seqnames = encode_blacklist_regions$chr,ranges = IRanges(start = encode_blacklist_regions$start+1,end = encode_blacklist_regions$end))
    gr_encode_blacklist <- GenomicRanges::reduce(gr_encode_blacklist)
    #用gem减去encode_blacklist
    gr_mappable <- GenomicRanges::setdiff(gr_gem_mappabilty, gr_encode_blacklist)}
  #2-1-7-2得到gr_vep_result(和2-2中的不一样，这里是1-based)
  if(1){
    seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
    #vep来自vcf，1-based
    start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
    end <- start
    gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
  hits_map_vepresult <- findOverlaps(gr_mappable,gr_vep_result)
  vep_result[,"Mappable"] <- FALSE
  vep_result[subjectHits(hits_map_vepresult),"Mappable"] <- TRUE
  
  #2-3-8添加 In_IR
  #2-3-8-1 IR本身就是GRanges对象了，不用获得
  #2-3-8-2得到gr_vep_result(和2-2中的不一样，这里是1-based)
  if(1){
    seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
    #vep来自vcf，1-based
    start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
    end <- start
    gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
  hits_IR_vep <- findOverlaps(IR,gr_vep_result)
  vep_result[,"IR"] <- FALSE
  vep_result[subjectHits(hits_IR_vep),"IR"] <- TRUE
  
  #2-3-9 添加靶基因expr的p和log2fc值(without cnv correction) #1h
  nc_rownames <- rownames(normalized_counts)
  nc_hugosymbol <- sapply(nc_rownames,function(x) strsplit(x,"\\|")[[1]][2])
  #result <- sapply(1:nrow(vep_result),function(i) p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol))
  cl <- makeCluster(40)
  registerDoParallel(cl)
  result <- foreach(i = 1:nrow(vep_result), .combine = c) %dopar% {
    p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol)
  }
  stopCluster(cl)
  vep_result$expr_p_log2fc <- result
  vep_result <- vep_result %>%
    separate(expr_p_log2fc,into = c("raw_expr_p","raw_expr_log2fc"),sep=":") 
    
    #2-3-10 添加靶基因expr的p_cnnormal和log2fc值(with cnv correction) #1h
    nc_rownames <- rownames(normalized_counts)
  nc_hugosymbol <- sapply(nc_rownames,function(x) strsplit(x,"\\|")[[1]][2])
  #result <- sapply(1:nrow(vep_result),function(i) p_expr(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol))
  cl <- makeCluster(40)
  registerDoParallel(cl)
  result <- foreach(i = 1:nrow(vep_result), .combine = c) %dopar% {
    library(rlang)
    p_expr_cnnormal(vep_result$Target_Gene[i], vep_result$Mutated_Samples[i],nc_rownames,nc_hugosymbol,cnv_1)
  }
  stopCluster(cl)
  vep_result$expr_p_log2fc_cnnormal <- result
  vep_result1 <- vep_result %>%
    separate(expr_p_log2fc_cnnormal,into = c("cnnormal_expr_p","cnnormal_expr_log2fc"),sep=":")
  #2-3-11 由于有一些基因在超过50%的基因中表达都为0，需要把这些基因算出来的expr_p去掉
  #得到这些基因在表达矩阵中对应的行名
  # 计算每行值为0的列的数量
  zero_count <- apply(normalized_counts == 0, 1, sum)
  # 计算总列数
  total_columns <- ncol(normalized_counts)
  # 筛选出小于50%列值为0的行的行名
  rownames  <- rownames(normalized_counts[zero_count < total_columns / 2, ])
  # 使用 strsplit 按照 "|" 分割
  split_vec <- strsplit(rownames, "\\|")
  # 扁平化为一个向量
  flattened_vec <- unlist(split_vec)
  #去掉vep_result1中的这些基因的expr
  vep_result1 <- vep_result1 %>%
    mutate(raw_expr_p=ifelse(Target_Gene %in% flattened_vec,raw_expr_p,NA)) %>%
    mutate(raw_expr_log2fc=ifelse(Target_Gene %in% flattened_vec,raw_expr_log2fc,NA)) %>%
    mutate(cnnormal_expr_p=ifelse(Target_Gene %in% flattened_vec,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_log2fc=ifelse(Target_Gene %in% flattened_vec,cnnormal_expr_log2fc,NA))
  #2-3-12 计算raw expr q值
  #raw_expr_p的NA可能是因为Func1有点问题，实际上是字符串
  vep_result1$raw_expr_p <- ifelse(!(vep_result1$raw_expr_p=="NA"),vep_result1$raw_expr_p,NA)
  vep_result1$raw_expr_p <- as.numeric(vep_result1$raw_expr_p)
  vep_result1$raw_expr_q<- p.adjust(vep_result1$raw_expr_p,method = "BH")
  #2-3-13 计算cnnormal expr q值
  vep_result1$cnnormal_expr_p <- as.numeric(vep_result1$cnnormal_expr_p)
  vep_result1$cnnormal_expr_q <- p.adjust(vep_result1$cnnormal_expr_p,method = "BH")
  enpm_output_10M <- vep_result1
  write.table(enpm_output_10M,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_output_10M.tsv",row.names = F,quote = F,sep = "\t")
}

#2-4 final_output
final_output <- enpm_output_tad[enpm_output_tad$Type !="Others",]
final_output$raw_expr_q <- p.adjust(final_output$raw_expr_p,method = "BH")
final_output$cnnormal_expr_q <- p.adjust(final_output$cnnormal_expr_p,method = "BH")
write.table(final_output,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_output.tsv",row.names = F,quote = F,sep = "\t")

#3添加得到各个突变的Mut Sig
#3-1 得到突变数据 
#导入vep注释文件
df_sig <- fread(enpmi_filepath)
df_sig$location <- paste(df_sig$V2,df_sig$V3,sep = ":")
df_sig <- df_sig[,-c(2,3)]
colnames(df_sig) <- c("Tumor_Sample_Barcode","ref","alt","location")
#去掉indel
df_sig <- df_sig %>%
  filter(nchar(ref) == 1 & nchar(alt) == 1) 
#设置一列key
df_sig <- df_sig %>%
  mutate(key = paste(Tumor_Sample_Barcode,location,sep = ":"))
#给all_singlemuts_sigs也加上key
all_singlemuts_sigs <- all_singlemuts_sigs %>%
  mutate(key = paste(Sample.Names,Chr,Pos,sep = ":"))
#得到Major_SBS
df_sig <- df_sig %>%
  left_join(all_singlemuts_sigs,by = "key") %>%
  select(Tumor_Sample_Barcode,ref,alt,location,Major_SBS)
#得到每个位置的突变数
df_sig <- df_sig %>%
  group_by(location) %>%
  mutate(mut.n = n()) %>%
  ungroup()
#给apobec和artifact添加注释
df_sig1 <- df_sig %>%
  mutate(Major_SBS = ifelse(Major_SBS %in% c("SBS2", "SBS13"), "APOBEC", Major_SBS)) %>%
  mutate(Major_SBS = ifelse(Major_SBS %in% paste("SBS",c(27, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 95),sep = ""), "Artifact", Major_SBS)) %>%
  mutate(Major_SBS = ifelse(Major_SBS %in% c("APOBEC","Artifact"), Major_SBS,"Others")) 
#统计一下是否有突变APOBEC+Artifact加起来>50%
#无
df_sig2 <- df_sig1 %>%
  # 创建一个新的列，标记 Major_SBS 是否是 "APOBEC" 或 "Artifact"
  mutate(is_APOBEC_or_Artifact = Major_SBS %in% c("APOBEC", "Artifact")) %>%
  # 按 location 分组
  group_by(location) %>%
  # 计算每个 location 中符合条件的比例
  summarize(
    count_total = n(),  # 计算每个 location 总的记录数
    count_APOBEC_or_Artifact = sum(is_APOBEC_or_Artifact)  # 计算符合条件的记录数
  ) %>%
  # 计算比例
  mutate(percentage = count_APOBEC_or_Artifact / count_total * 100) %>%
  # 筛选出比例大于等于 50% 的 location
  filter(percentage >= 50)

#4 绘图和表格筛选
#4-1 enpm_tad
enpm_tad_filepath <- file.path(folder_path,"enpm_output_tad.tsv")
df <- read.delim(enpm_tad_filepath)
df_x <- df[df$Type != "Others",]
df_x$raw_expr_q <- p.adjust(df_x$raw_expr_p,method = "BH")
df_x$cnnormal_expr_q <- p.adjust(df_x$cnnormal_expr_p,method = "BH")

#4-2 enpm_5M
enpm_5M_filepath <- file.path(folder_path,"enpm_output_5M.tsv")
df <- read.delim(enpm_5M_filepath)

#4-3 enpm_10M
enpm_10M_filepath <- file.path(folder_path,"enpm_output_10M.tsv")
df <- read.delim(enpm_10M_filepath)


if(0){
df <- vep_result1
#4-1 纵轴Location，横轴mut.n图 右边是expr_q
#寻找引起driver gene表达变化的点突变
#筛选条件1 类型为Promoter,Enhancer,Splice_site,UTR5,UTR3中的一种
#筛选条件2 靶基因为癌基因
df1  <- df %>%
  filter(Type %in% c("Enhancer","Promoter_CDS","SS_CDS","UTR5","UTR3")) %>%
  filter(Target_Gene %in% all_cancer_gene) %>%
  arrange(desc(mut.n))
write.table(df1,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_result_EPSU_cancergene",row.names = F,quote = F,sep = "\t")
#结论一:没有找到如TERT启动子上一样能明显引起基因表达变化的点突变!

#4-2 纵轴Location，横轴Samples图
#寻找引起RNA表达变化的点突变
#筛选条件1 类型为Promoter_LncRNA,SS_LncRNA
#筛选条件2 Target_gene q<0.1
df1 <- df %>%
  filter(Type %in% c("Promoter_LncRNA","SS_LncRNA")) %>%
  arrange(desc(mut.n))
draw_picture_q(df_sorted)
write.table(df1,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_result_epilncRNA",row.names = F,quote = F,sep = "\t")
#结论二：也没有找到能引起LncRNA表达明显变化的高频点突变


#4-3 纵轴Location，横轴Samples图
#寻找非编码RNA上的高频点突变
#筛选条件2 Target_gene是非编码RNA
df1 <- df %>%
  filter(Type %in% c("LncRNA","miRNA","SmallRNA")) %>%
  arrange(desc(mut.n))
write.table(df1,"/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_result_ncRNA",row.names = F,quote = F,sep = "\t")
#结论三：ncRNA上有高频的点突变，但意义不清楚

#4-4 纵轴Location，横轴Samples图
#寻找能引起基因表达改变的且不在注释元件上的高频点突变
#筛选条件1 类型为Others，即不在我找到的任何注释元件中
#筛选条件2 Target_gene expr q<0.1
df1 <- df %>%
  filter(Type %in% c("Others")) %>%
  filter(cnnormal_expr_q < 0.1) %>%
  arrange(desc(mut.n))
#无

#4-5 纵轴Location，横轴Samples图
#寻找和OS或RFS有关系的非编码点突变
df1 <- df %>%
  filter(os_p < 0.05 | rfs_p < 0.05) %>%
  arrange(desc(mut.n))
#有很多，但没想好怎么弄

#4-6 纵轴Location,横轴Samples,填充SBS类型
df <- df_sig1 %>%
  filter(mut.n > 15 & location!="chr7:55191822") %>%
  arrange(desc(mut.n))
df$Major_SBS <- factor(df$Major_SBS,levels=c("Others","APOBEC","Artifact"),ordered = T)
df_x <- df[,c("location","mut.n")]
df_x <- unique(df_x)
df$location <- factor(df$location, levels = df_x$location[order(df_x$mut.n)])
folder_path <- "/home/data/t190513/1000_noncoding/nc_point_mutations/"
filename <- file.path(folder_path,"Sigs_of_highfreq_singlemuts.pdf")
pdf(filename,width = 8,height = 6)
ggplot(df,aes(y=location,fill = Major_SBS)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = c("Others"= "grey", "APOBEC"= "#FDC086","Artifact"="#FFFF99"))
dev.off()

  

}

#草稿
