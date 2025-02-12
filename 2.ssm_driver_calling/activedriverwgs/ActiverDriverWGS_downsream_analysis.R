###ActiveDriverWGS 下游分析

#===============================================================================
##1 导入包，设置环境变量和自写函数
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(qqman)
library(detectIR)
library(stringr)
library(fdrtool)
library(doParallel)
library(GenomicRanges)
ADW_results_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/results_hg38_986jiaoji.tsv"
IR_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/IR_hg38.rds"
ADW_input_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/986_jiaoji.adwmuts.input"
all_elements_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_elements_hg38_noctcf.bed" 
GEM_mappability_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/GRCh38_mappability_100mer.gem.bed"
Encode_blacklist_regions_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/encode_hg38.blacklist.bed"
OncoKB_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/oncokb_cancerGeneList_2024_10_24.tsv"
CGC_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/Cancer Gene Census v100.tsv"
PCAWG_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/737 cancer genes by pcawg.csv"
lung_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/20220203_lung_drivers.csv"
counts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/counts_1000_log2tpmplus1.txt"
config_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/wgs_samples_986.txt"
cnv_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_data_by_genes.txt"
apobec_muts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/apobec_muts.txt"
TAD_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/A549_TAD_withgenes_hg38.bed"

#Func1 输入id，靶基因，和ems数据框，得到突变样本和非突变样本的值和检验p值和log2 FC
p_expr <- function(element_id,target_gene,df,expr_rownames,expr_hugosymbol){
  mutated_samples <- df %>% 
    filter(id==element_id) %>%
    select(sample)
  mutated_samples <- unique(mutated_samples$sample)
  mutated_indices <- which(names(expr) %in% mutated_samples)
  d_f <- data.frame()
  #如果该基因是ENSG开头,与在rownames里即可
  if(grepl("^ENSG",target_gene)){
    d_f <- expr[grepl(target_gene,expr_rownames),]
  }else{ #如果是Hugo_Symbol与Hugo_Symbol比对，需完全匹配
    d_f <- expr[expr_hugosymbol==target_gene,]
  }
  #如果能匹配上
  if(nrow(d_f)==1 & sum(mutated_indices) != 0){
    # 将表达量分组
    #原本的mutated samples是以，相连接的一些值，弄成一个向量
    group_mutated <- as.numeric(d_f[1,mutated_indices])
    group_unmutated <- as.numeric(d_f[1,-mutated_indices])
    #进行检验
    wilcox_test_result <- wilcox.test(group_mutated, group_unmutated)
    p <- wilcox_test_result$p.value
    #结果写入result
    log2_fc <- median(as.numeric(group_mutated))-median(as.numeric(group_unmutated))
    result <- paste(p,log2_fc,sep = ":")
    return(result)
  }else{return(NA)
  }
}
#Func2 输入id，靶基因，和ems数据框，得到突变样本和非突变样本的值和检验p值和log2 FC,和Func1不同的是
#假设检验仅在cn normal的样本中进行(1-22 cn=2)
p_expr_cnnormal <- function(element_id,target_gene,df,expr_rownames,expr_hugosymbol){
  mutated_samples <- df %>% 
    filter(id==element_id) %>%
    select(sample)
  mutated_samples <- unique(mutated_samples$sample)
  cn_samples <- names(cnv_1)[cnv_1[target_gene, ] == 0]
  cn_mutated_samples <- intersect(cn_samples,mutated_samples)
  cn_unmutated_samples <- setdiff(cn_samples,mutated_samples)
  mutated_indices <- which(names(expr) %in% cn_mutated_samples)
  unmutated_indices <- which(names(expr) %in% cn_unmutated_samples)
  d_f <- data.frame()
  #如果该基因是ENSG开头,与在rownames里即可
  if(grepl("^ENSG",target_gene)){
    d_f <- expr[grepl(target_gene,expr_rownames),]
  }else{ #如果是Hugo_Symbol与Hugo_Symbol比对，需完全匹配
    d_f <- expr[expr_hugosymbol==target_gene,]
  }
  #如果能匹配上且samples不为0
  if(nrow(d_f)==1 & sum(mutated_indices) != 0){
    # 将表达量分组
    #原本的mutated samples是以，相连接的一些值，弄成一个向量
    group_mutated <- as.numeric(d_f[1,mutated_indices])
    group_unmutated <- as.numeric(d_f[1,unmutated_indices])
    #进行检验
    wilcox_test_result <- wilcox.test(group_mutated, group_unmutated)
    p <- wilcox_test_result$p.value
    #结果写入result
    log2_fc <- median(as.numeric(group_mutated))-median(as.numeric(group_unmutated))
    result <- paste(p,log2_fc,sep = ":")
    return(result)
  }else{return(NA)
  }
}
#Func-3 绘图函数 输入基因的hugo symbol，返回共突变图
if(0){
#暂时还没弄好 draw_comutation <- function(gene_hugo_name){
  noncoding_MAF <- df[df$mut.n>=10 & df$Target_Gene == gene_hugo_name & (! is.na(df$Target_Gene))
                      & (df$expr_qh <0.1 | df$expr_ql < 0.1),]
  noncoding_MAF <- noncoding_MAF[apply(noncoding_MAF, 1, function(x) !all(is.na(x))), ]
  #设置一个id
  noncoding_MAF$id <- paste(noncoding_MAF$Target_Gene,noncoding_MAF$Location,sep = "_")
  #取出需要的一列并展开,弄成长数据
  noncoding_MAF <- noncoding_MAF[,c("id","Mutated_Samples")]
  noncoding_MAF <- noncoding_MAF %>%
    separate_rows(Mutated_Samples, sep = ",") %>%  # 拆分样本列
    mutate(value = 1)
  colnames(noncoding_MAF)[2] <- "Tumor_Sample_Barcode"
  noncoding_MAF <- unique(noncoding_MAF)
  #第二步 提取该基因编码区的突变
  coding_MAF <- coding_maf[coding_maf$Hugo_Symbol == gene_hugo_name,]
  if(nrow(coding_MAF > 0)){
    coding_MAF <- coding_MAF[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
    coding_MAF$Hugo_Symbol <- paste(coding_MAF$Hugo_Symbol,"_coding",sep="")
    coding_MAF$value <- 1
    colnames(coding_MAF)[1] <- "id"
    coding_MAF <- unique(coding_MAF)
    #第三步 合并MAF,转换成宽格式并画图
    MAF <- rbind(noncoding_MAF,coding_MAF)
    MAF <- MAF %>%
      pivot_wider(names_from = Tumor_Sample_Barcode, values_from = value,values_fill = 0)
    MAF <- as.data.frame(MAF)
    rownames(MAF) <- MAF$id
    MAF <- MAF[,-c(1)]
    mat <- as.matrix(MAF)
    col_fun = circlize::colorRamp2(c(0,1), c("grey","brown"))
    file_name <- paste(gene_hugo_name,"_comutation.jpg")
    jpeg(file_name,units = "cm",width = 10,height=4,res = 300)
    p <- Heatmap(mat,
                 col = col_fun,
                 cluster_rows = F,
                 cluster_columns = F,
                 rect_gp = gpar(col = "white", lwd = 3),
                 row_names_side = "left", 
                 show_column_names = F,
                 show_heatmap_legend = F,
                 border = T,
                 width = unit(4, "cm"), #设置单元格宽度
                 height = nrow(mat)*unit(6, "mm"),#设置单元格长度
    )
    draw(p)
    dev.off()
  }}

#=============================================================================== 
##2 导入和预处理文件
adw_results <- fread(ADW_results_filepath,data.table = F)
adw_results <- adw_results %>% 
  select(id,pp_element,element_muts_obs,element_muts_exp,element_enriched,fdr_element) %>%
  filter(! is.na(element_muts_exp))

#adw_results说明
#id 元素id
#pp_element 该元素突变富集的p值
#element_muts_obs 一共多少患者在该元素有突变
#element_muts_exp 按照背景突变概率该元素应该有多少突变
#element_enriched 富集(TRUE)还是不富集(FALSE)
#fdr_element 该元素突变富集的FDR q值
apobec_muts <- read.delim(apobec_muts_filepath) #1-based
apobec_muts <- apobec_muts%>%
  mutate(Chr = paste0("chr", Chr))
IR <- readRDS(IR_filepath) #1-based 双闭
adw_input_muts <- fread(ADW_input_filepath) #1-based 双闭
all_elements <- fread(all_elements_filepath) #0-based 左开右闭
names(all_elements) <- c("chr","start","end","id")
GEM_mappability <- fread(GEM_mappability_filepath) #0-based 左开右闭
names(GEM_mappability) <- c("chr","start","end")
encode_blacklist_regions <- fread(Encode_blacklist_regions_filepath) #0-based 左开右闭
names(encode_blacklist_regions) <- c("chr","start","end")


#得到一个需要用到的数据框 ems(element_mut_sample)
#具体说明
#id 该元件的id
#chr_ele 元件的chr
#start_ele 元件的start 1-based 双闭
#end_ele 元件的end 1-based 双闭
#chr_mut 突变的chr
#start_mut 突变的start 1-based 双闭
#end_mut 突变的end 1-based 双闭
#ref_mut 突变的reference  
#alt_mut 突变的alteration  
#sample 样本id
#In_IR 布尔值，该突变是否位于一段回文序列里
#Mappable 布尔值，该突变是否位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#APOBEC 布尔值，该突变是否来自于APOBEC突变过程
if(1){
#得到相应的granges对象
#本身就是1-based 双闭 直接转换成granges对象即可
gr_adw_input_muts <- GRanges(seqnames = adw_input_muts$chr,ranges = IRanges(start = adw_input_muts$pos1,end = adw_input_muts$pos2))
#是0-based 左闭右开 start+1 end+1-1
gr_all_elements <- GRanges(seqnames = all_elements$chr,ranges = IRanges(start = all_elements$start+1,end = all_elements$end))
#是0-based 左闭右开 start+1 end+1-1
gr_gem_mappabilty <- GRanges(seqnames = GEM_mappability$chr,ranges = IRanges(start = GEM_mappability$start+1,end = GEM_mappability$end))
gr_gem_mappabilty <- GenomicRanges::reduce(gr_gem_mappabilty)
#是0-based 左闭右开 start+1 end+1-1
gr_encode_blacklist <- GRanges(seqnames = encode_blacklist_regions$chr,ranges = IRanges(start = encode_blacklist_regions$start+1,end = encode_blacklist_regions$end))
gr_encode_blacklist <- GenomicRanges::reduce(gr_encode_blacklist)
#用gem减去encode_blacklist
gr_mappable <- GenomicRanges::setdiff(gr_gem_mappabilty, gr_encode_blacklist)
#apobec_muts 1-based 直接转换即可
gr_apobec_muts <- GRanges(seqnames = apobec_muts$Chr,ranges = IRanges(start = apobec_muts$Pos,end = apobec_muts$Pos))

#和element做Overlap 前面的是query,后面的是subject
hits_ele_mut <- findOverlaps(gr_adw_input_muts,gr_all_elements)
hits_ele_mut <- as.data.frame(hits_ele_mut)
ems <- adw_input_muts[hits_ele_mut$queryHits]
ems[,c("chr_ele","start_ele","end_ele","id")] <- all_elements[hits_ele_mut$subjectHits,]
#去掉NA的行(因为不是每个元件上都有突变)
ems <- na.omit(ems)
names(ems)[1:6] <- c("chr_mut","start_mut","end_mut","ref_mut","alt_mut","sample")
#和IR做Overlap
gr_emsmuts <- GRanges(seqnames=ems$chr_mut,IRanges(start = ems$start_mut,end = ems$end_mut))
hits_ems_mut <- findOverlaps(gr_emsmuts,IR)
ems[,"In_IR"] <- FALSE
ems[queryHits(hits_ems_mut),"In_IR"] <- TRUE

#和mappable做Overlap
hits_mut_map <- findOverlaps(gr_emsmuts,gr_mappable)
ems[,"Mappable"] <- FALSE
ems[queryHits(hits_mut_map),"Mappable"] <- TRUE

#和apobec做Overlap
hits_mut_apobec <- findOverlaps(gr_emsmuts,gr_apobec_muts)
ems[,"APOBEC"] <- FALSE
ems[queryHits(hits_mut_apobec),"APOBEC"] <- TRUE
#最后调整下all_elements的几项(因为是bed)
ems$start_ele <- ems$start_ele+1
#调整下列的顺序
ems <- ems %>%
  select(id,chr_ele,start_ele,end_ele,chr_mut,start_mut,end_mut,ref_mut,alt_mut,sample,In_IR,Mappable,APOBEC)

}

OncoKB_driver <- fread(OncoKB_drivergenes_filepath)
CGC_driver <- fread(CGC_drivergenes_filepath)
PCAWG_driver <- fread(PCAWG_drivergenes_filepath)
names(PCAWG_driver) <- "Gene_Symbol"
lung_driver <- fread(lung_drivergenes_filepath)
#得到config文件
config <- read.delim(config_filepath,header = F)
config$V1 <- sapply(config$V1, function(x) {str_remove(x, "(LC).*")})
#得到counts文件(log2(TPM+1))
expr <- fread(counts_filepath,data.table = F)
#fread不支持行名处理一下
rownames(expr) <- expr$V1
expr$V1 <- NULL
#只留下config有的样本Tumor的列
mid <- paste0(config$V1,"_T")
pattern <- paste(mid,collapse = "|")
expr <- expr[,grepl(pattern,colnames(expr))]
#改一下列名
colnames(expr) <- sapply(colnames(expr),function(x) paste0(strsplit(x,"_")[[1]][3],strsplit(x,"_")[[1]][2]))

expr_rownames <- rownames(expr)
expr_hugosymbol <- sapply(expr_rownames,function(x) strsplit(x,"\\|")[[1]][2])

#拷贝数结果
cnv <- fread(cnv_filepath,data.table = F) #里面的值是log2(CN)-1
cnv_1 <- cnv
rownames(cnv_1) <- cnv_1$`Gene Symbol`
cnv_1[,1:3] <- NULL
#所有的癌基因
all_drivergenes_list <- union(OncoKB_driver$`Hugo Symbol`,CGC_driver$`Gene Symbol`)
all_drivergenes_list  <- union(all_drivergenes_list,PCAWG_driver$Gene_Symbol)             
all_drivergenes_list <- union(all_drivergenes_list,lung_driver$Gene_Symbol)

#所有的肺癌 癌基因
CGC_NSCLC_drivergene_list <- CGC_driver %>%
  filter(str_detect("NSCLC",`Tumour Types(Somatic)`)) %>%
  select(`Gene Symbol`)
lung_drivergenes_list <- union(lung_driver$Gene_Symbol,CGC_NSCLC_drivergene_list$`Gene Symbol`)

#===============================================================================
##3 质控
#3-1 qqplot
#跑图 要个几分钟
qq(adw_results$pp_element)
#计算膨胀系数 预期为1 越偏离1说明群体分层现象严重，容易有假阳性
p_value <- adw_results$pp_element
z <- qnorm(p_value/2) 
lambda  <- round(median(z^2, na.rm = TRUE) / 0.456, 3)
#lambda的值为0.422

##3-2 post-filtering
#3-2-1 至少3个患者有该突变 去掉了约1/4
adw_results <- adw_results %>%
  mutate(element_muts_obs_morethan2 = element_muts_obs >= 3)

#3-2-2 less than 50% mutations are located in palindromic DNA sequence
#得到超过50%突变在回文序列的元件id 去除了8000多个
id_IRratio_overhalf <- ems %>%
  group_by(id) %>%
  summarise(
    total_muts_n = n(),
    IR_muts_n = sum(In_IR == TRUE)
  ) %>%
  mutate(IRratio = IR_muts_n / total_muts_n) %>%
  filter(IRratio >= 0.5) %>%
  select(id) 
adw_results <- adw_results %>%
  mutate(IRratio_overhalf = id %in% id_IRratio_overhalf$id)
#3-2-3 more than 50% of mutations are located in mappable genomic regions (CRGalignability(100mer), 
#DAC blacklisted regions)
#得到超过50%突变在mappable区域的元件id 去掉了3000多个
id_Mapratio_overhalf <- ems %>%
  group_by(id) %>%
  summarise(
    total_muts_n = n(),
    mappable_muts_n = sum(Mappable == TRUE)
  ) %>%
  mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
  filter(Mapratio >= 0.5) %>%
  select(id) 
adw_results <- adw_results %>%
  mutate(Mapratio_overhalf = id %in% id_Mapratio_overhalf$id)
#3-2-4 less than 50% of mutations attributed to APOBEC mutation signatures
#计算每个元件突变是否小于50%来源于APOBEC
id_APOBEC_lesshalf <- ems %>%
  group_by(id) %>%
  summarise(
    total_muts_n = n(),
    apobec_muts_n = sum(APOBEC == TRUE)
  ) %>%
  mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
  filter(apobecratio < 0.5) %>%
  select(id) 
adw_results <- adw_results %>%
  mutate(Apobecratio_lesshalf = id %in% id_APOBEC_lesshalf$id)

#3-2-5 在筛除之后把筛掉的元件的p设为1，重新计算q值
adw_results <- adw_results %>%
  mutate(Post_filter_p = case_when(
    IRratio_overhalf == FALSE & Mapratio_overhalf == TRUE & element_muts_obs >2 & Apobecratio_lesshalf == TRUE ~ pp_element ,
    TRUE ~ 1))
fdr <- fdrtool(adw_results$Post_filter_p,statistic="pvalue")
adw_results$Post_filter_q <- fdr$qval

#3-3 计算sensitivity和precision
#利用找到的candidate coding gene占已知的lung driver gene比例反映sensitivity
#利用找到的candidate coding gene中已知driver gene(任何癌症)的比例反映precision
candidate_coding_element <- adw_results %>%
  filter(Post_filter_q < 0.1 & str_detect(id,"cds"))
ccg_list <- sapply(candidate_coding_element$id,function(x) strsplit(x,"::")[[1]][3])
sen <- sum(lung_drivergenes_list %in% ccg_list) / length(lung_drivergenes_list)
#sensitivity = 0.185
pre <- sum(ccg_list %in% all_drivergenes_list) / length(ccg_list)
#precision = 0.1

#4 计算调控元件对靶基因的expr effect，得到最终表格
#表格组成

#添加Type
if(1){
adw_results <- adw_results %>%
  mutate(Type = case_when(
    str_detect(id, "gc46_pc.cds") ~ "CDS",
    str_detect(id, "gc46_pc.ss_cds") ~ "SS_CDS",
    str_detect(id, "gc46_pc.5utr") ~ "UTR5",
    str_detect(id, "gc46_pc.3utr") ~ "UTR3",
    str_detect(id, "gc46.promoter") ~ "Promoter_CDS",
    str_detect(id, "lncrna_exons") ~ "LncRNA",
    str_detect(id, "ss_lncrna") ~ "SS_LncRNA",
    str_detect(id, "lncrna.promoter") ~ "Promoter_LncRNA",
    str_detect(id, "smallrna") ~ "SmallRNA",
    str_detect(id, "mirna") ~ "miRNA",
    str_detect(id, "enhancers") ~ "Enhancer",
    TRUE ~ "Others"  # 其他情况赋值为 Others
  ))  }
#添加Target_gene
if(1){
adw_results$Target_gene <- NA
rows <- str_detect(adw_results$id,"enhancers") 
#非enhancer元件
adw_results[-which(rows),"Target_gene"] <- sapply(adw_results$id[-which(rows)],function(x) strsplit(x,"::")[[1]][3])
#enhancer
adw_results[rows,"Target_gene"] <- sapply(adw_results$id[rows],function(x) strsplit(x,"::")[[1]][4])
}
#添加raw_expr_p及raw_expr_log2fc(median)
if(1){
results <- rep(NA,nrow(adw_results))
cl <- makeCluster(20)
# 注册并行后端
registerDoParallel(cl)
# 并行执行 for 循环
results <- foreach(i = 1:nrow(adw_results), .combine = c) %dopar% {
  library(dplyr)
  p_expr(adw_results$id[i], adw_results$Target_gene[i], ems, expr_rownames, expr_hugosymbol)
}
# 停止并行计算
stopCluster(cl)
adw_results$p_expr <- results
adw_results <- adw_results %>%
  separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
}

#添加cnnormal_expr_p及cnnormal_expr_log2fc
if(1){
  results <- rep(NA,nrow(adw_results))
  cl <- makeCluster(20)
  # 注册并行后端
  registerDoParallel(cl)
  # 并行执行 for 循环
  results <- foreach(i = 1:nrow(adw_results), .combine = c) %dopar% {
    library(dplyr)
    p_expr_cnnormal(adw_results$id[i], adw_results$Target_gene[i], ems, expr_rownames, expr_hugosymbol)
  }
  # 停止并行计算
  stopCluster(cl)
  adw_results$p_expr_cnnormal <- results
  adw_results <- adw_results %>%
    separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
}

#去除<50%肿瘤中表达基因的expr_p和log2FC
if(1){
# 计算每行值为0的列的数量
zero_count <- apply(expr == 0, 1, sum)
# 计算总列数
total_columns <- ncol(expr)
# 筛选出小于50%列值为0的行的行名
rownames  <- rownames(expr[zero_count < total_columns / 2, ])
# 使用 strsplit 按照 "|" 分割
split_vec <- strsplit(rownames, "\\|")
# 扁平化为一个向量
flattened_vec <- unlist(split_vec)
#去掉adw_results中的这些基因的expr_p和log2FC
adw_results <- adw_results %>%
  mutate(raw_expr_p=ifelse(Target_gene %in% flattened_vec,raw_expr_p,NA)) %>%
  mutate(raw_expr_log2FC=ifelse(Target_gene %in% flattened_vec,raw_expr_log2FC,NA)) %>%
  mutate(cnnormal_expr_p=ifelse(Target_gene %in% flattened_vec,cnnormal_expr_p,NA)) %>%
  mutate(cnnormal_expr_log2FC=ifelse(Target_gene %in% flattened_vec,cnnormal_expr_log2FC,NA))
}

#计算raw_expr_q和cnnormal_expr_q
adw_results$raw_expr_p <- as.numeric(adw_results$raw_expr_p)
fdr <- fdrtool(adw_results$raw_expr_p[! is.na(adw_results$raw_expr_p)],statistic="pvalue")
adw_results$raw_expr_q[! is.na(adw_results$raw_expr_p)] <- fdr$qval
#计算cnnormal expr q值
adw_results$cnnormal_expr_p <- as.numeric(adw_results$cnnormal_expr_p)
fdr <- fdrtool(abs(adw_results$cnnormal_expr_p[! is.na(adw_results$cnnormal_expr_p)]),statistic="pvalue")
adw_results$cnnormal_expr_q[! is.na(adw_results$cnnormal_expr_p)] <- fdr$qval

write.table(adw_results,file ="/home/data/t190513/1000_noncoding/activedriverwgs/adw_results_jiaoji_2025_0119.txt",row.names = F,quote = F,sep = "\t")

#==================================================================================================



