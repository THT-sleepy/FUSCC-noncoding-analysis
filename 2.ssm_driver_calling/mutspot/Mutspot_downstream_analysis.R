#MutSpot下游分析
#由于tnseq 突变的结果一直跑不好，先用取交集的突变进行分析

#1 导入需要的包和自写函数以及文件路径
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(GenomicRanges)
library(doParallel)
library(qqman)
folder_path <- "/home/data/t190513/1000_noncoding/mutspot/MutSpot/"
snv_hotspots_merged_filepath <- file.path(folder_path,"986_jiaoji_result/snv_hotspots_merged.tsv")
indel_hotspots_merged_filepath <- file.path(folder_path,"986_jiaoji_result/indel_hotspots_merged.tsv")
all_elements_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_elements_hg38_noctcf.bed" 
tad_withgenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/A549_TAD_withgenes_hg38.bed"
counts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/counts_1000_log2tpmplus1.txt"
config_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/wgs_samples_986.txt"
snv_input_filepath <- file.path(folder_path,"986_jiaoji.mutspot_snv.input")
indel_input_filepath <- file.path(folder_path,"986_jiaoji.mutspot_indel.input")
GEM_mappability_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/GRCh38_mappability_100mer.gem.bed"
Encode_blacklist_regions_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/encode_hg38.blacklist.bed"
IR_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/IR_hg38.rds"
cnv_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_data_by_genes.txt"
OncoKB_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/oncokb_cancerGeneList_2024_10_24.tsv"
CGC_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/Cancer Gene Census v100.tsv"
PCAWG_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/737 cancer genes by pcawg.csv"
lung_drivergenes_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/20220203_lung_drivers.csv"
all_gene_bed_filepath <- file.path(folder_path,"all_gene_hg38.bed4")
apobec_muts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/apobec_muts.txt"
#Func1 输入id，靶基因，和hss/hls数据框，以及expr矩阵的行名及对应的hugosymbol(这两项是每次循环会用到的变量
#，单独设置就不用每次循环都再设置，速度更快)得到突变样本和非突变样本的值和检验p值和log2 FC
p_expr <- function(hotspot,target_gene,df,expr_rownames,expr_hugosymbol){
  id <- hotspot
  mutated_samples <- df %>% 
    filter(hotspot==id) %>%
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
#Func2 输入同Func1 和Func1不同的是假设检验仅在cn normal的样本中进行(1-22 cn=2,n=820)
p_expr_cnnormal <- function(hotspot,target_gene,df,expr_rownames,expr_hugosymbol){
  id <- hotspot
  mutated_samples <- df %>% 
    filter(hotspot==id) %>%
    select(sample)
  mutated_samples <- unique(mutated_samples$sample)
  cn_samples <- names(cnv)[cnv[target_gene, ] == 0]
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
#2 导入和预处理文件
if(1){
#2-1 snv hotspots merged
#id
#chrom 
#start(1-based)
#end(1-based)
#pval p value
#length length of hotspot
#p.bg Mean background mutation probability
#k number of mutated samples
#fdr fdr矫正后的q值
snv_hotspots_merged <- fread(snv_hotspots_merged_filepath,data.table = F)
rownames(snv_hotspots_merged) <- snv_hotspots_merged$V1
snv_hotspots_merged$V1 <- NULL
#取出q<0.1的
snv_hotspots_merged <- snv_hotspots_merged %>%
  mutate(id=rownames(snv_hotspots_merged))%>%
  filter(fdr<0.1) 
  

#2-2 indel hotspots merged
#id
#chrom 
#start(1-based)
#end(1-based)
#pval p value
#length length of hotspot
#p.bg Mean background mutation probability
#k number of mutated samples
#fdr fdr矫正后的q值
indel_hotspots_merged <- fread(indel_hotspots_merged_filepath,data.table = F)
rownames(indel_hotspots_merged) <- indel_hotspots_merged$V1
indel_hotspots_merged$V1 <- NULL
#取出q<0.1的
indel_hotspots_merged <- indel_hotspots_merged %>%
  mutate(id=rownames(indel_hotspots_merged))%>%
  filter(fdr<0.1) 
  

#2-3 所有注释元件的坐标文件(0-based 左闭右开)
all_elements <- fread(all_elements_filepath) #0-based 左闭右开


#2-4 tad的坐标文件(0-based 左闭右开)
tad <- fread(tad_withgenes_filepath)

#2-5 expr文件 log2(TPM+1)以及config文件
config <- read.delim(config_filepath,header = F)
config$V1 <- sapply(config$V1, function(x) {str_remove(x, "(LC).*")})
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

#2-6 input_muts(1-based，左闭右闭)
input_snv <- fread(snv_input_filepath)
names(input_snv) <- c("chr","start","end","ref","alt","sample")
input_indel <- fread(indel_input_filepath)
names(input_indel) <- c("chr","start","end","ref","alt","sample")
#input_muts <- rbind(input_snv,input_indel)

#2-7 GEM mappability(0-based,左闭右开)
GEM_mappability <- fread(GEM_mappability_filepath)
names(GEM_mappability) <- c("chr","start","end")

#2-8 Encode blacklist (0-based,左闭右开)
encode_blacklist_regions <- fread(Encode_blacklist_regions_filepath) 
names(encode_blacklist_regions) <- c("chr","start","end")

#2-9 APOBEC muts(1-based)
apobec_muts <- read.delim(apobec_muts_filepath) #1-based
apobec_muts <- apobec_muts%>%
  mutate(Chr = paste0("chr", Chr))

#2-10 IR Inverse Repeats
IR <- readRDS(IR_filepath) #1-based 双闭

#2-11 拷贝数变异结果
cnv <- fread(cnv_filepath,data.table = F) #里面的值是log2(CN)-1
cnv <- cnv
rownames(cnv) <- cnv$`Gene Symbol`
cnv[,1:3] <- NULL

#2-12 癌症相关基因
OncoKB_driver <- fread(OncoKB_drivergenes_filepath)
CGC_driver <- fread(CGC_drivergenes_filepath)
PCAWG_driver <- fread(PCAWG_drivergenes_filepath)
names(PCAWG_driver) <- "Gene_Symbol"
lung_driver <- fread(lung_drivergenes_filepath)
#所有的癌基因
all_drivergenes_list <- union(OncoKB_driver$`Hugo Symbol`,CGC_driver$`Gene Symbol`)
all_drivergenes_list  <- union(all_drivergenes_list,PCAWG_driver$Gene_Symbol)             
all_drivergenes_list <- union(all_drivergenes_list,lung_driver$Gene_Symbol)
#所有的肺癌 癌基因
CGC_NSCLC_drivergene_list <- CGC_driver %>%
  filter(str_detect("NSCLC",`Tumour Types(Somatic)`)) %>%
  select(`Gene Symbol`)
lung_drivergenes_list <- union(lung_driver$Gene_Symbol,CGC_NSCLC_drivergene_list$`Gene Symbol`)

#2-13 all_gene 坐标区间(0-based 左闭右开)
all_gene_interval <- fread(all_gene_bed_filepath)
names(all_gene_interval) <- c("chr","start","end","gene")
}

#3 
#说明
#3-1 #得到hotspot_snv_sample(hss)数据框
if(0){
#由于hotspot是以区间为单位进行描述的，丢失了每个突变相应的信息，
#而在计算mutation-expression association(需要突变对应的样本)，以及
#进行质控，如统计是否>50%的突变位于apobec区间时，都需要用到这些
#信息，所以和ADW_downsream_analysis一样，需要生成一个以单个突变为单位
#进行描述的表格
}
#id 该hotspot的id
#chr_hotspot hotspot的chr
#start_hotspot hotspot的start 1-based 双闭
#end_hotspot hotspot的end 1-based 双闭
#chr_mut 突变的chr
#start_mut 突变的start 1-based 双闭
#end_mut 突变的end 1-based 双闭
#ref_mut 突变的reference  
#alt_mut 突变的alteration  
#sample 突变对应的样本id
#In_IR 布尔值，该突变是否位于一段回文序列里
#Mappable 布尔值，该突变是否位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#APOBEC，该突变是否很可能来源于APOBEC介导的突变过程
if(1){
  #得到相应的granges对象
  #本身就是1-based 双闭 直接转换成granges对象即可
  gr_input_snv <- GRanges(seqnames = input_snv$chr,ranges = IRanges(start = input_snv$start,end = input_snv$end))
  #是1-based 双闭，左不变，右不变
  gr_hotspot <- GRanges(seqnames = snv_hotspots_merged$chr,ranges = IRanges(start = snv_hotspots_merged$start,end = snv_hotspots_merged$end))
  #是0-based 左闭右开 start+1 end+1-1
  gr_gem_mappabilty <- GRanges(seqnames = GEM_mappability$chr,ranges = IRanges(start = GEM_mappability$start+1,end = GEM_mappability$end))
  gr_gem_mappabilty <- GenomicRanges::reduce(gr_gem_mappabilty)
  #是0-based 左闭右开 start+1 end+1-1
  gr_encode_blacklist <- GRanges(seqnames = encode_blacklist_regions$chr,ranges = IRanges(start = encode_blacklist_regions$start+1,end = encode_blacklist_regions$end))
  gr_encode_blacklist <- GenomicRanges::reduce(gr_encode_blacklist)
  #用gem减去encode_blacklist
  gr_mappable <- GenomicRanges::setdiff(gr_gem_mappabilty, gr_encode_blacklist)
  #apobec_muts 1-based 直接转换即可
  #gr_apobec_muts <- GRanges(seqnames = apobec_muts$Chr,ranges = IRanges(start = apobec_muts$Pos,end = apobec_muts$Pos))
  #和snv hotspot(sh)做Overlap 前面的是query,后面的是subject
  hits_mut_sh <- findOverlaps(gr_input_snv,gr_hotspot)
  hits_mut_sh <- as.data.frame(hits_mut_sh)
  hss <- input_snv[hits_mut_sh$queryHits]
  hss[,c("chr_hotspot","start_hotspot","end_hotspot","hotspot")] <- snv_hotspots_merged[hits_mut_sh$subjectHits,c("chrom","start","end","id")]
  names(hss)[1:6] <- c("chr_mut","start_mut","end_mut","ref_mut","alt_mut","sample")
  #和IR做Overlap
  gr_hssmuts <- GRanges(seqnames=hss$chr_mut,IRanges(start = hss$start_mut,end = hss$end_mut))
  hits_hss_mut <- findOverlaps(gr_hssmuts,IR)
  hss[,"In_IR"] <- FALSE
  hss[queryHits(hits_hss_mut),"In_IR"] <- TRUE
  
  #和mappable做Overlap
  hits_mut_map <- findOverlaps(gr_hssmuts,gr_mappable)
  hss[,"Mappable"] <- FALSE
  hss[queryHits(hits_mut_map),"Mappable"] <- TRUE
  
  #和apobec做Overlap
  if(1){
  gr_apobec_muts <- GRanges(seqnames = apobec_muts$Chr,ranges = IRanges(start=apobec_muts$Pos,end = apobec_muts$Pos))
  hits_mut_apobec <- findOverlaps(gr_hssmuts,gr_apobec_muts)
  hss[,"APOBEC"] <- FALSE
  hss[queryHits(hits_mut_apobec),"APOBEC"] <- TRUE}
  #调整下列的顺序
  hss <- hss %>%
    select(hotspot,chr_hotspot,start_hotspot,end_hotspot,chr_mut,start_mut,end_mut,ref_mut,alt_mut,sample,In_IR,Mappable,APOBEC)
  
}

#3-2 #得到hotspot_indel_sample(his)数据框
#id 该hotspot的id
#chr_hotspot hotspot的chr
#start_hotspot hotspot的start 1-based 双闭
#end_hotspot hotspot的end 1-based 双闭
#chr_mut 突变的chr
#start_mut 突变的start 1-based 双闭
#end_mut 突变的end 1-based 双闭
#ref_mut 突变的reference  
#alt_mut 突变的alteration  
#sample 突变对应的样本id
#In_IR 布尔值，该突变是否位于一段回文序列里
#Mappable 布尔值，该突变是否位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#APOBEC，该突变是否很可能来源于APOBEC介导的突变过程
if(1){
  #得到相应的granges对象
  #本身就是1-based 双闭 直接转换成granges对象即可
  gr_input_indel <- GRanges(seqnames = input_indel$chr,ranges = IRanges(start = input_indel$start,end = input_indel$end))
  #是1-based 双闭，左不变，右不变
  gr_hotspot <- GRanges(seqnames = indel_hotspots_merged$chr,ranges = IRanges(start = indel_hotspots_merged$start,end = indel_hotspots_merged$end))
  #是0-based 左闭右开 start+1 end+1-1
  gr_gem_mappabilty <- GRanges(seqnames = GEM_mappability$chr,ranges = IRanges(start = GEM_mappability$start+1,end = GEM_mappability$end))
  gr_gem_mappabilty <- GenomicRanges::reduce(gr_gem_mappabilty)
  #是0-based 左闭右开 start+1 end+1-1
  gr_encode_blacklist <- GRanges(seqnames = encode_blacklist_regions$chr,ranges = IRanges(start = encode_blacklist_regions$start+1,end = encode_blacklist_regions$end))
  gr_encode_blacklist <- GenomicRanges::reduce(gr_encode_blacklist)
  #用gem减去encode_blacklist
  gr_mappable <- GenomicRanges::setdiff(gr_gem_mappabilty, gr_encode_blacklist)
  #apobec_muts 1-based 直接转换即可
  #gr_apobec_muts <- GRanges(seqnames = apobec_muts$Chr,ranges = IRanges(start = apobec_muts$Pos,end = apobec_muts$Pos))
  #和indel hotspot(sh)做Overlap 前面的是query,后面的是subject
  hits_mut_sh <- findOverlaps(gr_input_indel,gr_hotspot)
  hits_mut_sh <- as.data.frame(hits_mut_sh)
  his <- input_indel[hits_mut_sh$queryHits]
  his[,c("chr_hotspot","start_hotspot","end_hotspot","hotspot")] <- indel_hotspots_merged[hits_mut_sh$subjectHits,c("chrom","start","end","id")]
  names(his)[1:6] <- c("chr_mut","start_mut","end_mut","ref_mut","alt_mut","sample")
  #和IR做Overlap
  gr_hismuts <- GRanges(seqnames=his$chr_mut,IRanges(start = his$start_mut,end = his$end_mut))
  hits_his_mut <- findOverlaps(gr_hismuts,IR)
  his[,"In_IR"] <- FALSE
  his[queryHits(hits_his_mut),"In_IR"] <- TRUE
  
  #和mappable做Overlap
  hits_mut_map <- findOverlaps(gr_hismuts,gr_mappable)
  his[,"Mappable"] <- FALSE
  his[queryHits(hits_mut_map),"Mappable"] <- TRUE
  
  #和apobec做Overlap
  if(1){
    hits_mut_apobec <- findOverlaps(gr_hismuts,gr_apobec_muts)
    his[,"APOBEC"] <- FALSE
    his[queryHits(hits_mut_apobec),"APOBEC"] <- TRUE}
  #调整下列的顺序
  his <- his %>%
    select(hotspot,chr_hotspot,start_hotspot,end_hotspot,chr_mut,start_mut,end_mut,ref_mut,alt_mut,sample,In_IR,Mappable,APOBEC) 
  
}


#4 计算mutation-expression association
#4-1 snv,within注释元件取对应基因，others取相同TAD基因
if(1){
#最终得到数据框snv_tad
#hotspot 该hotspot对应的id
#chrom 
#start(1-based)
#end(1-based)
#pval p value
#length length of hotspot
#p.bg Mean background mutation probability
#k number of mutated samples
#fdr fdr矫正后的q值
#element within的元件(hotspot within元件)
#type 
#target_gene 假定的靶基因,within注释元件即为对应基因，others为相同TAD的某个基因
#raw_expr_p 突变组相较于非突变组的wilconxin p
#raw_expr_log2FC
#cnnormal_expr_p 突变组相较于非突变组的wilconxin p,仅在拷贝数正常的样本中进行比较
#cnnormal_expr_log2FC
#raw_expr_allgene_q 在所有基因(>50%样本中表达)中进行fdr 矫正
#raw_expr_driveronly_q 仅在drivergene(>50%样本中表达)中进行fdr 矫正
#cnnormal_expr_allgene_q
#cnnormal_expr_driveronly_q
#IRratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于一段回文序列里
#Mapratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
#Apobecratio_lesshalf 布尔值，该hotspot是否小于百分之50突变来自于APOBEC突变过程
#4-1-1 添加overlapping element,type和注释元件的target_gene
if(1){
#0-based ，左闭右开,左+1,右不变
gr_ele <- GRanges(seqnames = all_elements$V1,ranges = IRanges(start = all_elements$V2+1,end = all_elements$V3))
#1-based，左闭右闭，左闭右闭
gr_snvhotspots <- GRanges(seqnames = snv_hotspots_merged$chrom,ranges = IRanges(start = snv_hotspots_merged$start,end = snv_hotspots_merged$end))
#做overlap 前面的是query,后面的是subject,query within subject
hits_ele_snvhotspots <- findOverlaps(gr_snvhotspots,gr_ele,type="within")
#由于一个hotspot可对应多个element，需要整合一下(1个hotspot的element弄到一行)
overlap <- data.frame(element=all_elements$V4[subjectHits(hits_ele_snvhotspots)],hotspot=rownames(snv_hotspots_merged)[queryHits(hits_ele_snvhotspots)])
overlap <- overlap %>%
  #添加type
  mutate(type = case_when(
    str_detect(element, "gc46_pc.cds") ~ "CDS",
    str_detect(element, "gc46_pc.ss_cds") ~ "SS_CDS",
    str_detect(element, "gc46_pc.5utr") ~ "UTR5",
    str_detect(element, "gc46_pc.3utr") ~ "UTR3",
    str_detect(element, "gc46.promoter") ~ "Promoter_CDS",
    str_detect(element, "lncrna_exons") ~ "LncRNA",
    str_detect(element, "ss_lncrna") ~ "SS_LncRNA",
    str_detect(element, "lncrna.promoter") ~ "Promoter_LncRNA",
    str_detect(element, "smallrna") ~ "SmallRNA",
    str_detect(element, "mirna") ~ "miRNA",
    str_detect(element, "enhancers") ~ "Enhancer",
    TRUE ~ "Others"  # 其他情况赋值为 Others
  )) %>%
  #添加target_gene
  mutate(target_gene = ifelse(!(type %in% c("Others","Enhancer")),sapply(overlap$element,function(x) strsplit(x,"::")[[1]][3]),
                          ifelse(type == "Enhancer",sapply(overlap$element,function(x) strsplit(x,"::")[[1]][4]),NA))) %>%
  #将每个hotspot对应的多个element整合到一起
  mutate(key1=paste(element,target_gene,type,sep = ":::")) %>%
  select(hotspot,key1) %>%
  group_by(hotspot) %>%
  summarise(key2=paste(key1,collapse = ";"))
snv_tad <- snv_hotspots_merged %>%
  mutate(hotspot=rownames(snv_hotspots_merged)) %>%
  left_join(overlap,by="hotspot") %>%
  separate_rows(key2,sep = ";") %>%
  separate(key2,into = c("element","target_gene","type"),sep = ":::") %>%
  mutate(type=ifelse(is.na(type),"Others",type))
}
#4-1-2 添加注释元件外的target_gene
if(1){
#0-based 左闭右开，左+1，右不变
gr_tad <- GRanges(seqnames = tad$chr,ranges = IRanges(start = tad$start+1,end = tad$end))
#1-based 左闭右闭，左不变，右不变
gr_snv_tad <- GRanges(seqnames = snv_tad$chrom,ranges = IRanges(start = snv_tad$start,end = snv_tad$end))
hits_tad_snvtad <- findOverlaps(query = gr_tad,subject = gr_snv_tad)
#经检查，每个hotspot最多仅对应一个TAD
#添加others的target_gene
snv_tad$others_target[subjectHits(hits_tad_snvtad)] <- tad$genes_within[queryHits(hits_tad_snvtad)]
snv_tad <- snv_tad %>%
  mutate(target_gene=ifelse(is.na(target_gene),others_target,target_gene)) %>%
  separate_rows(target_gene,sep = ";")
#去掉所有target_gene在小于50%样本中表达的行
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
snv_tad <- snv_tad %>%
  filter(target_gene %in% flattened_vec)
}
}
#4-1-3 添加raw_expr_p,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_log2FC
#添加raw_expr_p,raw_expr_log2FC
if(1){
  results <- rep(NA,nrow(snv_tad))
  cl <- makeCluster(40)
  # 注册并行后端
  registerDoParallel(cl)
  # 并行执行 for 循环
  results <- foreach(i = 1:nrow(snv_tad), .combine = c) %dopar% {
    library(dplyr)
    p_expr(snv_tad$hotspot[i], snv_tad$target_gene[i], hss, expr_rownames, expr_hugosymbol)
  }
  # 停止并行计算
  stopCluster(cl)
  snv_tad$p_expr <- results
  snv_tad <- snv_tad %>%
    separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
}
#添加cnnormal_expr_p,cnnormal_expr_log2FC
if(1){
  results <- rep(NA,nrow(snv_tad))
  cl <- makeCluster(40)
  # 注册并行后端
  registerDoParallel(cl)
  # 并行执行 for 循环
  results <- foreach(i = 1:nrow(snv_tad), .combine = c) %dopar% {
    library(dplyr)
    p_expr_cnnormal(snv_tad$hotspot[i], snv_tad$target_gene[i], hss, expr_rownames, expr_hugosymbol)
  }
  # 停止并行计算
  stopCluster(cl)
  snv_tad$p_expr_cnnormal <- results
  snv_tad <- snv_tad %>%
    separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
}
#添加#raw_expr_allgene_q,raw_expr_driveronly_q,cnnormal_expr_allgene_q,cnnormal_expr_driveronly_q
snv_tad <- snv_tad %>%
  mutate(raw_expr_allgene_q=p.adjust(raw_expr_p,method = "BH")) %>%
  mutate(raw_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,raw_expr_p,NA)) %>%
  mutate(raw_expr_driveronly_q=p.adjust(raw_expr_driveronly_p,method = "BH")) %>%
  mutate(cnnormal_expr_allgene_q=p.adjust(cnnormal_expr_p,method = "BH")) %>%
  mutate(cnnormal_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,cnnormal_expr_p,NA)) %>%
  mutate(cnnormal_expr_driveronly_q=p.adjust(cnnormal_expr_driveronly_p,method = "BH")) %>%
  select(-raw_expr_driveronly_p,-cnnormal_expr_driveronly_p,-others_target)
#4-1-4添加IRratio_overhalf,Mapratio_overhalf,Apobecratio_lesshalf
if(1){
#IRratio_overhalf
hotspot_IRratio_overhalf <- hss %>%
  group_by(hotspot) %>%
  summarise(
    total_muts_n = n(),
    IR_muts_n = sum(In_IR == TRUE)
  ) %>%
  mutate(IRratio = IR_muts_n / total_muts_n) %>%
  filter(IRratio >= 0.5) %>%
  select(hotspot) 
snv_tad <- snv_tad %>%
  mutate(IRratio_overhalf = hotspot %in% hotspot_IRratio_overhalf$hotspot)
#Mapratio_overhalf
hotspot_Mapratio_overhalf <- hss %>%
  group_by(hotspot) %>%
  summarise(
    total_muts_n = n(),
    mappable_muts_n = sum(Mappable == TRUE)
  ) %>%
  mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
  filter(Mapratio >= 0.5) %>%
  select(hotspot) 
snv_tad <- snv_tad %>%
  mutate(Mapratio_overhalf = hotspot %in% hotspot_Mapratio_overhalf$hotspot)
#Apobecratio_lesshalf
hotspot_APOBEC_lesshalf <- hss %>%
  group_by(hotspot) %>%
  summarise(
    total_muts_n = n(),
    apobec_muts_n = sum(APOBEC == TRUE)
  ) %>%
  mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
  filter(apobecratio < 0.5) %>%
  select(hotspot) 
snv_tad <- snv_tad %>%
  mutate(Apobecratio_lesshalf = hotspot %in% hotspot_APOBEC_lesshalf$hotspot)}
filename <- file.path(folder_path,"snv_tad.txt")
write.table(snv_tad,file =filename,row.names = F,quote = F,sep = "\t")
}

#4-2 snv,within注释元件取对应基因，others取距离<5M基因
if(1){
  #最终得到数据框snv_5M
  #hotspot 该hotspot对应的id
  #chrom 
  #start(1-based)
  #end(1-based)
  #pval p value
  #length length of hotspot
  #p.bg Mean background mutation probability
  #k number of mutated samples
  #fdr fdr矫正后的q值
  #element within的元件(hotspot within元件)
  #type 
  #target_gene 假定的靶基因,within注释元件即为对应基因，others为距离5M内的某个基因
  #raw_expr_p 突变组相较于非突变组的wilconxin p
  #raw_expr_log2FC
  #cnnormal_expr_p 突变组相较于非突变组的wilconxin p,仅在拷贝数正常的样本中进行比较
  #cnnormal_expr_log2FC
  #raw_expr_allgene_q 在所有基因(>50%样本中表达)中进行fdr 矫正
  #raw_expr_driveronly_q 仅在drivergene(>50%样本中表达)中进行fdr 矫正
  #cnnormal_expr_allgene_q
  #cnnormal_expr_driveronly_q
  #IRratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于一段回文序列里
  #Mapratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
  #Apobecratio_lesshalf 布尔值，该hotspot是否小于百分之50突变来自于APOBEC突变过程
  
  #4-2-1 添加overlapping element,type和overlap element的target_gene
  if(1){
    #0-based ，左闭右开,左+1,右不变
    gr_ele <- GRanges(seqnames = all_elements$V1,ranges = IRanges(start = all_elements$V2+1,end = all_elements$V3))
    #1-based，左闭右闭，左闭右闭
    gr_snvhotspots <- GRanges(seqnames = snv_hotspots_merged$chrom,ranges = IRanges(start = snv_hotspots_merged$start,end = snv_hotspots_merged$end))
    #做overlap 前面的是query,后面的是subject,query within subject
    hits_ele_snvhotspots <- findOverlaps(gr_snvhotspots,gr_ele,type="within")
    #由于一个hotspot可对应多个element，需要整合一下(1个hotspot的element弄到一行)
    overlap <- data.frame(element=all_elements$V4[subjectHits(hits_ele_snvhotspots)],hotspot=rownames(snv_hotspots_merged)[queryHits(hits_ele_snvhotspots)])
    overlap <- overlap %>%
      #添加type
      mutate(type = case_when(
        str_detect(element, "gc46_pc.cds") ~ "CDS",
        str_detect(element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(element, "lncrna_exons") ~ "LncRNA",
        str_detect(element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(element, "smallrna") ~ "SmallRNA",
        str_detect(element, "mirna") ~ "miRNA",
        str_detect(element, "enhancers") ~ "Enhancer",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      )) %>%
      #添加target_gene
      mutate(target_gene = ifelse(!(type %in% c("Others","Enhancer")),sapply(overlap$element,function(x) strsplit(x,"::")[[1]][3]),
                                  ifelse(type == "Enhancer",sapply(overlap$element,function(x) strsplit(x,"::")[[1]][4]),NA))) %>%
      #将每个hotspot对应的多个element整合到一起
      mutate(key1=paste(element,target_gene,type,sep = ":::")) %>%
      select(hotspot,key1) %>%
      group_by(hotspot) %>%
      summarise(key2=paste(key1,collapse = ";"))
    snv_5M <- snv_hotspots_merged %>%
      mutate(hotspot=rownames(snv_hotspots_merged)) %>%
      left_join(overlap,by="hotspot") %>%
      separate_rows(key2,sep = ";") %>%
      separate(key2,into = c("element","target_gene","type"),sep = ":::") %>%
      mutate(type=ifelse(is.na(type),"Others",type))
  }
  #4-2-2 添加others的target_gene
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_gene <- GRanges(seqnames = all_gene_interval$chr,ranges = IRanges(start = all_gene_interval$start+1,end = all_gene_interval$end))
    #1-based 左闭右闭，左不变，右不变
    gr_snv_5M <- GRanges(seqnames = snv_5M$chrom,ranges = IRanges(start = snv_5M$start,end = snv_5M$end))
    hits_gene_snv5M <- findOverlaps(query = gr_gene,subject = gr_snv_5M,maxgap = 5000000)
    #添加非pro/utr的target_gene
    #每个非pro/utr的hotspot对应多个target_gene,先把它们整理成每个hotspot一行
    hotspot_target <- data.frame(hotspot=snv_5M$hotspot[subjectHits(hits_gene_snv5M)],others_target=all_gene_interval$gene[queryHits(hits_gene_snv5M)])
    hotspot_target <- hotspot_target %>%
      group_by(hotspot) %>%
      summarise(others_targets=paste(others_target,collapse = ";"))
    snv_5M <- snv_5M %>%
      left_join(hotspot_target,by="hotspot") %>%
      mutate(target_gene=ifelse(is.na(target_gene),others_targets,target_gene)) %>%
      separate_rows(target_gene,sep = ";")
    #只保留target_gene在大于50%样本中表达的行
    if(1){
      # 计算每行值为0的列的数量
      zero_count <- apply(expr == 0, 1, sum)
      # 计算总列数
      total_columns <- ncol(expr)
      # 筛选出大于50%列值为0的行的行名
      rownames  <- rownames(expr[zero_count < total_columns / 2, ])
      # 使用 strsplit 按照 "|" 分割
      split_vec <- strsplit(rownames, "\\|")
      # 扁平化为一个向量
      flattened_vec <- unlist(split_vec)
      snv_5M <- snv_5M %>%
        filter(target_gene %in% flattened_vec)
    }
  }
  #4-2-3 添加raw_expr_p,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_log2FC
  #添加raw_expr_p,raw_expr_log2FC
  if(1){
    results <- rep(NA,nrow(snv_5M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(snv_5M), .combine = c) %dopar% {
      library(dplyr)
      p_expr(snv_5M$hotspot[i], snv_5M$target_gene[i], hss, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    snv_5M$p_expr <- results
    snv_5M <- snv_5M %>%
      separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
  }
  #添加cnnormal_expr_p,cnnormal_expr_log2FC
  if(1){
    results <- rep(NA,nrow(snv_5M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(snv_5M), .combine = c) %dopar% {
      library(dplyr)
      p_expr_cnnormal(snv_5M$hotspot[i], snv_5M$target_gene[i], hss, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    snv_5M$p_expr_cnnormal <- results
    snv_5M <- snv_5M %>%
      separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
  }
  #添加#raw_expr_allgene_q,raw_expr_driveronly_q,cnnormal_expr_allgene_q,cnnormal_expr_driveronly_q
  snv_5M <- snv_5M %>%
    mutate(raw_expr_allgene_q=p.adjust(raw_expr_p,method = "BH")) %>%
    mutate(raw_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,raw_expr_p,NA)) %>%
    mutate(raw_expr_driveronly_q=p.adjust(raw_expr_driveronly_p,method = "BH")) %>%
    mutate(cnnormal_expr_allgene_q=p.adjust(cnnormal_expr_p,method = "BH")) %>%
    mutate(cnnormal_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_driveronly_q=p.adjust(cnnormal_expr_driveronly_p,method = "BH")) %>%
    select(-raw_expr_driveronly_p,-cnnormal_expr_driveronly_p,-others_targets)
  #4-2-4添加IRratio_overhalf,Mapratio_overhalf,Apobecratio_lesshalf
  if(1){
    #IRratio_overhalf
    hotspot_IRratio_overhalf <- hss %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        IR_muts_n = sum(In_IR == TRUE)
      ) %>%
      mutate(IRratio = IR_muts_n / total_muts_n) %>%
      filter(IRratio >= 0.5) %>%
      select(hotspot) 
    snv_5M <- snv_5M %>%
      mutate(IRratio_overhalf = hotspot %in% hotspot_IRratio_overhalf$hotspot)
    #Mapratio_overhalf
    hotspot_Mapratio_overhalf <- hss %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        mappable_muts_n = sum(Mappable == TRUE)
      ) %>%
      mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
      filter(Mapratio >= 0.5) %>%
      select(hotspot) 
    snv_5M <- snv_5M %>%
      mutate(Mapratio_overhalf = hotspot %in% hotspot_Mapratio_overhalf$hotspot)
    #Apobecratio_lesshalf
    hotspot_APOBEC_lesshalf <- hss %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        apobec_muts_n = sum(APOBEC == TRUE)
      ) %>%
      mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
      filter(apobecratio < 0.5) %>%
      select(hotspot) 
    snv_5M <- snv_5M %>%
      mutate(Apobecratio_lesshalf = hotspot %in% hotspot_APOBEC_lesshalf$hotspot)}
  filename <- file.path(folder_path,"snv_5M.txt")
  write.table(snv_5M,file =filename,row.names = F,quote = F,sep = "\t")
}

#4-3 snv,within注释元件取对应基因，others取距离<10M基因
if(1){
  #最终得到数据框snv_10M
  #hotspot 该hotspot对应的id
  #chrom 
  #start(1-based)
  #end(1-based)
  #pval p value
  #length length of hotspot
  #p.bg Mean background mutation probability
  #k number of mutated samples
  #fdr fdr矫正后的q值
  #element within的元件(hotspot within元件)
  #type 
  #target_gene 假定的靶基因,within注释元件即为对应基因，others为距离10M的某个基因
  #raw_expr_p 突变组相较于非突变组的wilconxin p
  #raw_expr_log2FC
  #cnnormal_expr_p 突变组相较于非突变组的wilconxin p,仅在拷贝数正常的样本中进行比较
  #cnnormal_expr_log2FC
  #raw_expr_allgene_q 在所有基因(>50%样本中表达)中进行fdr 矫正
  #raw_expr_driveronly_q 仅在drivergene(>50%样本中表达)中进行fdr 矫正
  #cnnormal_expr_allgene_q
  #cnnormal_expr_driveronly_q
  #IRratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于一段回文序列里
  #Mapratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
  #Apobecratio_lesshalf 布尔值，该hotspot是否小于百分之50突变来自于APOBEC突变过程
  
  #4-3-1 添加overlapping element,type和overlap element的target_gene
  if(1){
    #0-based ，左闭右开,左+1,右不变
    gr_ele <- GRanges(seqnames = all_elements$V1,ranges = IRanges(start = all_elements$V2+1,end = all_elements$V3))
    #1-based，左闭右闭，左闭右闭
    gr_snvhotspots <- GRanges(seqnames = snv_hotspots_merged$chrom,ranges = IRanges(start = snv_hotspots_merged$start,end = snv_hotspots_merged$end))
    #做overlap 前面的是query,后面的是subject,query within subject
    hits_ele_snvhotspots <- findOverlaps(gr_snvhotspots,gr_ele,type="within")
    #由于一个hotspot可对应多个element，需要整合一下(1个hotspot的element弄到一行)
    overlap <- data.frame(element=all_elements$V4[subjectHits(hits_ele_snvhotspots)],hotspot=rownames(snv_hotspots_merged)[queryHits(hits_ele_snvhotspots)])
    overlap <- overlap %>%
      #添加type
      mutate(type = case_when(
        str_detect(element, "gc46_pc.cds") ~ "CDS",
        str_detect(element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(element, "lncrna_exons") ~ "LncRNA",
        str_detect(element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(element, "smallrna") ~ "SmallRNA",
        str_detect(element, "mirna") ~ "miRNA",
        str_detect(element, "enhancers") ~ "Enhancer",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      )) %>%
      #添加target_gene
      mutate(target_gene = ifelse(!(type %in% c("Others","Enhancer")),sapply(overlap$element,function(x) strsplit(x,"::")[[1]][3]),
                                  ifelse(type == "Enhancer",sapply(overlap$element,function(x) strsplit(x,"::")[[1]][4]),NA))) %>%
      #将每个hotspot对应的多个element整合到一起
      mutate(key1=paste(element,target_gene,type,sep = ":::")) %>%
      select(hotspot,key1) %>%
      group_by(hotspot) %>%
      summarise(key2=paste(key1,collapse = ";"))
    snv_10M <- snv_hotspots_merged %>%
      mutate(hotspot=rownames(snv_hotspots_merged)) %>%
      left_join(overlap,by="hotspot") %>%
      separate_rows(key2,sep = ";") %>%
      separate(key2,into = c("element","target_gene","type"),sep = ":::") %>%
      mutate(type=ifelse(is.na(type),"Others",type))
  }
  #4-3-2 添加others的target_gene
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_gene <- GRanges(seqnames = all_gene_interval$chr,ranges = IRanges(start = all_gene_interval$start+1,end = all_gene_interval$end))
    #1-based 左闭右闭，左不变，右不变
    gr_snv_10M <- GRanges(seqnames = snv_10M$chrom,ranges = IRanges(start = snv_10M$start,end = snv_10M$end))
    hits_gene_snv10M <- findOverlaps(query = gr_gene,subject = gr_snv_10M,maxgap = 10000000)
    #添加非pro/utr的target_gene
    #每个非pro/utr的hotspot对应多个target_gene,先把它们整理成每个hotspot一行
    hotspot_target <- data.frame(hotspot=snv_10M$hotspot[subjectHits(hits_gene_snv10M)],others_target=all_gene_interval$gene[queryHits(hits_gene_snv10M)])
    hotspot_target <- hotspot_target %>%
      group_by(hotspot) %>%
      summarise(others_targets=paste(others_target,collapse = ";"))
    snv_10M <- snv_10M %>%
      left_join(hotspot_target,by="hotspot") %>%
      mutate(target_gene=ifelse(is.na(target_gene),others_targets,target_gene)) %>%
      separate_rows(target_gene,sep = ";")
    #只保留target_gene在大于50%样本中表达的行
    if(1){
      # 计算每行值为0的列的数量
      zero_count <- apply(expr == 0, 1, sum)
      # 计算总列数
      total_columns <- ncol(expr)
      # 筛选出大于50%列值为0的行的行名
      rownames  <- rownames(expr[zero_count < total_columns / 2, ])
      # 使用 strsplit 按照 "|" 分割
      split_vec <- strsplit(rownames, "\\|")
      # 扁平化为一个向量
      flattened_vec <- unlist(split_vec)
      snv_10M <- snv_10M %>%
        filter(target_gene %in% flattened_vec)
    }
  }
  #4-3-3 添加raw_expr_p,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_log2FC
  #添加raw_expr_p,raw_expr_log2FC
  if(1){
    results <- rep(NA,nrow(snv_10M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(snv_10M), .combine = c) %dopar% {
      library(dplyr)
      p_expr(snv_10M$hotspot[i], snv_10M$target_gene[i], hss, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    snv_10M$p_expr <- results
    snv_10M <- snv_10M %>%
      separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
  }
  #添加cnnormal_expr_p,cnnormal_expr_log2FC
  if(1){
    results <- rep(NA,nrow(snv_10M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(snv_10M), .combine = c) %dopar% {
      library(dplyr)
      p_expr_cnnormal(snv_10M$hotspot[i], snv_10M$target_gene[i], hss, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    snv_10M$p_expr_cnnormal <- results
    snv_10M <- snv_10M %>%
      separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
  }
  #添加#raw_expr_allgene_q,raw_expr_driveronly_q,cnnormal_expr_allgene_q,cnnormal_expr_driveronly_q
  snv_10M <- snv_10M %>%
    mutate(raw_expr_allgene_q=p.adjust(raw_expr_p,method = "BH")) %>%
    mutate(raw_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,raw_expr_p,NA)) %>%
    mutate(raw_expr_driveronly_q=p.adjust(raw_expr_driveronly_p,method = "BH")) %>%
    mutate(cnnormal_expr_allgene_q=p.adjust(cnnormal_expr_p,method = "BH")) %>%
    mutate(cnnormal_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_driveronly_q=p.adjust(cnnormal_expr_driveronly_p,method = "BH")) %>%
    select(-raw_expr_driveronly_p,-cnnormal_expr_driveronly_p,-others_targets)
  #4-3-4添加IRratio_overhalf,Mapratio_overhalf,Apobecratio_lesshalf
  if(1){
    #IRratio_overhalf
    hotspot_IRratio_overhalf <- hss %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        IR_muts_n = sum(In_IR == TRUE)
      ) %>%
      mutate(IRratio = IR_muts_n / total_muts_n) %>%
      filter(IRratio >= 0.5) %>%
      select(hotspot) 
    snv_10M <- snv_10M %>%
      mutate(IRratio_overhalf = hotspot %in% hotspot_IRratio_overhalf$hotspot)
    #Mapratio_overhalf
    hotspot_Mapratio_overhalf <- hss %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        mappable_muts_n = sum(Mappable == TRUE)
      ) %>%
      mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
      filter(Mapratio >= 0.5) %>%
      select(hotspot) 
    snv_10M <- snv_10M %>%
      mutate(Mapratio_overhalf = hotspot %in% hotspot_Mapratio_overhalf$hotspot)
    #Apobecratio_lesshalf
    hotspot_APOBEC_lesshalf <- hss %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        apobec_muts_n = sum(APOBEC == TRUE)
      ) %>%
      mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
      filter(apobecratio < 0.5) %>%
      select(hotspot) 
    snv_10M <- snv_10M %>%
      mutate(Apobecratio_lesshalf = hotspot %in% hotspot_APOBEC_lesshalf$hotspot)}
  filename <- file.path(folder_path,"snv_10M.txt")
  write.table(snv_10M,file =filename,row.names = F,quote = F,sep = "\t")
}

#4-4 indel,within注释元件取对应基因，others取相同TAD基因
if(1){
  #最终得到数据框indel_tad
  #hotspot 该hotspot对应的id
  #chrom 
  #start(1-based)
  #end(1-based)
  #pval p value
  #length length of hotspot
  #p.bg Mean background mutation probability
  #k number of mutated samples
  #fdr fdr矫正后的q值
  #element within的元件(hotspot within元件)
  #type 
  #target_gene 假定的靶基因,within注释元件即为对应基因，others为相同TAD的某个基因
  #raw_expr_p 突变组相较于非突变组的wilconxin p
  #raw_expr_log2FC
  #cnnormal_expr_p 突变组相较于非突变组的wilconxin p,仅在拷贝数正常的样本中进行比较
  #cnnormal_expr_log2FC
  #raw_expr_allgene_q 在所有基因(>50%样本中表达)中进行fdr 矫正
  #raw_expr_driveronly_q 仅在drivergene(>50%样本中表达)中进行fdr 矫正
  #cnnormal_expr_allgene_q
  #cnnormal_expr_driveronly_q
  #IRratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于一段回文序列里
  #Mapratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
  #Apobecratio_lesshalf 布尔值，该hotspot是否小于百分之50突变来自于APOBEC突变过程
  
  #4-4-1 添加overlapping element,type和注释元件的target_gene
  if(1){
    #0-based ，左闭右开,左+1,右不变
    gr_ele <- GRanges(seqnames = all_elements$V1,ranges = IRanges(start = all_elements$V2+1,end = all_elements$V3))
    #1-based，左闭右闭，左闭右闭
    gr_indelhotspots <- GRanges(seqnames = indel_hotspots_merged$chrom,ranges = IRanges(start = indel_hotspots_merged$start,end = indel_hotspots_merged$end))
    #做overlap 前面的是query,后面的是subject,query within subject
    hits_ele_indelhotspots <- findOverlaps(gr_indelhotspots,gr_ele,type="within")
    #由于一个hotspot可对应多个element，需要整合一下(1个hotspot的element弄到一行)
    overlap <- data.frame(element=all_elements$V4[subjectHits(hits_ele_indelhotspots)],hotspot=rownames(indel_hotspots_merged)[queryHits(hits_ele_indelhotspots)])
    overlap <- overlap %>%
      #添加type
      mutate(type = case_when(
        str_detect(element, "gc46_pc.cds") ~ "CDS",
        str_detect(element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(element, "lncrna_exons") ~ "LncRNA",
        str_detect(element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(element, "smallrna") ~ "SmallRNA",
        str_detect(element, "mirna") ~ "miRNA",
        str_detect(element, "enhancers") ~ "Enhancer",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      )) %>%
      #添加target_gene
      mutate(target_gene = ifelse(!(type %in% c("Others","Enhancer")),sapply(overlap$element,function(x) strsplit(x,"::")[[1]][3]),
                                  ifelse(type == "Enhancer",sapply(overlap$element,function(x) strsplit(x,"::")[[1]][4]),NA))) %>%
      #将每个hotspot对应的多个element整合到一起
      mutate(key1=paste(element,target_gene,type,sep = ":::")) %>%
      select(hotspot,key1) %>%
      group_by(hotspot) %>%
      summarise(key2=paste(key1,collapse = ";"))
    indel_tad <- indel_hotspots_merged %>%
      mutate(hotspot=rownames(indel_hotspots_merged)) %>%
      left_join(overlap,by="hotspot") %>%
      separate_rows(key2,sep = ";") %>%
      separate(key2,into = c("element","target_gene","type"),sep = ":::") %>%
      mutate(type=ifelse(is.na(type),"Others",type))
  }
  #4-4-2 添加非pro/utr的target_gene
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_tad <- GRanges(seqnames = tad$chr,ranges = IRanges(start = tad$start+1,end = tad$end))
    #1-based 左闭右闭，左不变，右不变
    gr_indel_tad <- GRanges(seqnames = indel_tad$chrom,ranges = IRanges(start = indel_tad$start,end = indel_tad$end))
    hits_tad_indeltad <- findOverlaps(query = gr_tad,subject = gr_indel_tad)
    #经检查，每个hotspot最多仅对应一个TAD
    #添加非pro/utr的target_gene
    indel_tad$others_target[subjectHits(hits_tad_indeltad)] <- tad$genes_within[queryHits(hits_tad_indeltad)]
    indel_tad <- indel_tad %>%
      mutate(target_gene=ifelse(is.na(target_gene),others_target,target_gene)) %>%
      separate_rows(target_gene,sep = ";")
    #去掉所有target_gene在小于50%样本中表达的行
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
      indel_tad <- indel_tad %>%
        filter(target_gene %in% flattened_vec)
    }
  }
  #4-4-3 添加raw_expr_p,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_log2FC
  #添加raw_expr_p,raw_expr_log2FC
  if(1){
    results <- rep(NA,nrow(indel_tad))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(indel_tad), .combine = c) %dopar% {
      library(dplyr)
      p_expr(indel_tad$hotspot[i], indel_tad$target_gene[i], his, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    indel_tad$p_expr <- results
    indel_tad <- indel_tad %>%
      separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
  }
  #添加cnnormal_expr_p,cnnormal_expr_log2FC
  if(1){
    results <- rep(NA,nrow(indel_tad))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(indel_tad), .combine = c) %dopar% {
      library(dplyr)
      p_expr_cnnormal(indel_tad$hotspot[i], indel_tad$target_gene[i], his, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    indel_tad$p_expr_cnnormal <- results
    indel_tad <- indel_tad %>%
      separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
  }
  #添加#raw_expr_allgene_q,raw_expr_driveronly_q,cnnormal_expr_allgene_q,cnnormal_expr_driveronly_q
  indel_tad <- indel_tad %>%
    mutate(raw_expr_allgene_q=p.adjust(raw_expr_p,method = "BH")) %>%
    mutate(raw_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,raw_expr_p,NA)) %>%
    mutate(raw_expr_driveronly_q=p.adjust(raw_expr_driveronly_p,method = "BH")) %>%
    mutate(cnnormal_expr_allgene_q=p.adjust(cnnormal_expr_p,method = "BH")) %>%
    mutate(cnnormal_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_driveronly_q=p.adjust(cnnormal_expr_driveronly_p,method = "BH")) %>%
    select(-raw_expr_driveronly_p,-cnnormal_expr_driveronly_p,-others_target)
  #4-4-4添加IRratio_overhalf,Mapratio_overhalf,Apobecratio_lesshalf
  if(1){
    #IRratio_overhalf
    hotspot_IRratio_overhalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        IR_muts_n = sum(In_IR == TRUE)
      ) %>%
      mutate(IRratio = IR_muts_n / total_muts_n) %>%
      filter(IRratio >= 0.5) %>%
      select(hotspot) 
    indel_tad <- indel_tad %>%
      mutate(IRratio_overhalf = hotspot %in% hotspot_IRratio_overhalf$hotspot)
    #Mapratio_overhalf
    hotspot_Mapratio_overhalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        mappable_muts_n = sum(Mappable == TRUE)
      ) %>%
      mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
      filter(Mapratio >= 0.5) %>%
      select(hotspot) 
    indel_tad <- indel_tad %>%
      mutate(Mapratio_overhalf = hotspot %in% hotspot_Mapratio_overhalf$hotspot)
    #Apobecratio_lesshalf
    hotspot_APOBEC_lesshalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        apobec_muts_n = sum(APOBEC == TRUE)
      ) %>%
      mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
      filter(apobecratio < 0.5) %>%
      select(hotspot) 
    indel_tad <- indel_tad %>%
      mutate(Apobecratio_lesshalf = hotspot %in% hotspot_APOBEC_lesshalf$hotspot)}
  filename <- file.path(folder_path,"indel_tad.txt")
  write.table(indel_tad,file =filename,row.names = F,quote = F,sep = "\t")
}

#4-5 indel,within注释元件取对应基因，others取距离<5M基因
if(1){
  #最终得到数据框indel_5M
  #hotspot 该hotspot对应的id
  #chrom 
  #start(1-based)
  #end(1-based)
  #pval p value
  #length length of hotspot
  #p.bg Mean background mutation probability
  #k number of mutated samples
  #fdr fdr矫正后的q值
  #element within的元件(hotspot within元件)
  #type 
  #target_gene 假定的靶基因,within注释元件即为对应基因，others为距离5M内的某个基因
  #raw_expr_p 突变组相较于非突变组的wilconxin p
  #raw_expr_log2FC
  #cnnormal_expr_p 突变组相较于非突变组的wilconxin p,仅在拷贝数正常的样本中进行比较
  #cnnormal_expr_log2FC
  #raw_expr_allgene_q 在所有基因(>50%样本中表达)中进行fdr 矫正
  #raw_expr_driveronly_q 仅在drivergene(>50%样本中表达)中进行fdr 矫正
  #cnnormal_expr_allgene_q
  #cnnormal_expr_driveronly_q
  #IRratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于一段回文序列里
  #Mapratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
  #Apobecratio_lesshalf 布尔值，该hotspot是否小于百分之50突变来自于APOBEC突变过程
  
  #4-5-1 添加overlapping element,type和注释元件的target_gene
  if(1){
    #0-based ，左闭右开,左+1,右不变
    gr_ele <- GRanges(seqnames = all_elements$V1,ranges = IRanges(start = all_elements$V2+1,end = all_elements$V3))
    #1-based，左闭右闭，左闭右闭
    gr_indelhotspots <- GRanges(seqnames = indel_hotspots_merged$chrom,ranges = IRanges(start = indel_hotspots_merged$start,end = indel_hotspots_merged$end))
    #做overlap 前面的是query,后面的是subject,query within subject
    hits_ele_indelhotspots <- findOverlaps(gr_indelhotspots,gr_ele,type="within")
    #由于一个hotspot可对应多个element，需要整合一下(1个hotspot的element弄到一行)
    overlap <- data.frame(element=all_elements$V4[subjectHits(hits_ele_indelhotspots)],hotspot=rownames(indel_hotspots_merged)[queryHits(hits_ele_indelhotspots)])
    overlap <- overlap %>%
      #添加type
      mutate(type = case_when(
        str_detect(element, "gc46_pc.cds") ~ "CDS",
        str_detect(element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(element, "lncrna_exons") ~ "LncRNA",
        str_detect(element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(element, "smallrna") ~ "SmallRNA",
        str_detect(element, "mirna") ~ "miRNA",
        str_detect(element, "enhancers") ~ "Enhancer",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      )) %>%
      #添加target_gene
      mutate(target_gene = ifelse(!(type %in% c("Others","Enhancer")),sapply(overlap$element,function(x) strsplit(x,"::")[[1]][3]),
                                  ifelse(type == "Enhancer",sapply(overlap$element,function(x) strsplit(x,"::")[[1]][4]),NA))) %>%
      #将每个hotspot对应的多个element整合到一起
      mutate(key1=paste(element,target_gene,type,sep = ":::")) %>%
      select(hotspot,key1) %>%
      group_by(hotspot) %>%
      summarise(key2=paste(key1,collapse = ";"))
    indel_5M <- indel_hotspots_merged %>%
      mutate(hotspot=rownames(indel_hotspots_merged)) %>%
      left_join(overlap,by="hotspot") %>%
      separate_rows(key2,sep = ";") %>%
      separate(key2,into = c("element","target_gene","type"),sep = ":::") %>%
      mutate(type=ifelse(is.na(type),"Others",type))
  }
  #4-5-2 添加非pro/utr的target_gene
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_gene <- GRanges(seqnames = all_gene_interval$chr,ranges = IRanges(start = all_gene_interval$start+1,end = all_gene_interval$end))
    #1-based 左闭右闭，左不变，右不变
    gr_indel_5M <- GRanges(seqnames = indel_5M$chrom,ranges = IRanges(start = indel_5M$start,end = indel_5M$end))
    hits_gene_indel5M <- findOverlaps(query = gr_gene,subject = gr_indel_5M,maxgap = 5000000)
    #添加非pro/utr的target_gene
    #每个非pro/utr的hotspot对应多个target_gene,先把它们整理成每个hotspot一行
    hotspot_target <- data.frame(hotspot=indel_5M$hotspot[subjectHits(hits_gene_indel5M)],others_target=all_gene_interval$gene[queryHits(hits_gene_indel5M)])
    hotspot_target <- hotspot_target %>%
      group_by(hotspot) %>%
      summarise(others_targets=paste(others_target,collapse = ";"))
    indel_5M <- indel_5M %>%
      left_join(hotspot_target,by="hotspot") %>%
      mutate(target_gene=ifelse(is.na(target_gene),others_targets,target_gene)) %>%
      separate_rows(target_gene,sep = ";")
    #只保留target_gene在大于50%样本中表达的行
    if(1){
      # 计算每行值为0的列的数量
      zero_count <- apply(expr == 0, 1, sum)
      # 计算总列数
      total_columns <- ncol(expr)
      # 筛选出大于50%列值为0的行的行名
      rownames  <- rownames(expr[zero_count < total_columns / 2, ])
      # 使用 strsplit 按照 "|" 分割
      split_vec <- strsplit(rownames, "\\|")
      # 扁平化为一个向量
      flattened_vec <- unlist(split_vec)
      indel_5M <- indel_5M %>%
        filter(target_gene %in% flattened_vec)
    }
  }
  #4-5-3 添加raw_expr_p,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_log2FC
  #添加raw_expr_p,raw_expr_log2FC
  if(1){
    results <- rep(NA,nrow(indel_5M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(indel_5M), .combine = c) %dopar% {
      library(dplyr)
      p_expr(indel_5M$hotspot[i], indel_5M$target_gene[i], his, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    indel_5M$p_expr <- results
    indel_5M <- indel_5M %>%
      separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
  }
  #添加cnnormal_expr_p,cnnormal_expr_log2FC
  if(1){
    results <- rep(NA,nrow(indel_5M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(indel_5M), .combine = c) %dopar% {
      library(dplyr)
      p_expr_cnnormal(indel_5M$hotspot[i], indel_5M$target_gene[i], his, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    indel_5M$p_expr_cnnormal <- results
    indel_5M <- indel_5M %>%
      separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
  }
  #添加#raw_expr_allgene_q,raw_expr_driveronly_q,cnnormal_expr_allgene_q,cnnormal_expr_driveronly_q
  indel_5M <- indel_5M %>%
    mutate(raw_expr_allgene_q=p.adjust(raw_expr_p,method = "BH")) %>%
    mutate(raw_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,raw_expr_p,NA)) %>%
    mutate(raw_expr_driveronly_q=p.adjust(raw_expr_driveronly_p,method = "BH")) %>%
    mutate(cnnormal_expr_allgene_q=p.adjust(cnnormal_expr_p,method = "BH")) %>%
    mutate(cnnormal_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_driveronly_q=p.adjust(cnnormal_expr_driveronly_p,method = "BH")) %>%
    select(-raw_expr_driveronly_p,-cnnormal_expr_driveronly_p,-others_targets)
  #4-5-4添加IRratio_overhalf,Mapratio_overhalf,Apobecratio_lesshalf
  if(1){
    #IRratio_overhalf
    hotspot_IRratio_overhalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        IR_muts_n = sum(In_IR == TRUE)
      ) %>%
      mutate(IRratio = IR_muts_n / total_muts_n) %>%
      filter(IRratio >= 0.5) %>%
      select(hotspot) 
    indel_5M <- indel_5M %>%
      mutate(IRratio_overhalf = hotspot %in% hotspot_IRratio_overhalf$hotspot)
    #Mapratio_overhalf
    hotspot_Mapratio_overhalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        mappable_muts_n = sum(Mappable == TRUE)
      ) %>%
      mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
      filter(Mapratio >= 0.5) %>%
      select(hotspot) 
    indel_5M <- indel_5M %>%
      mutate(Mapratio_overhalf = hotspot %in% hotspot_Mapratio_overhalf$hotspot)
    #Apobecratio_lesshalf
    hotspot_APOBEC_lesshalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        apobec_muts_n = sum(APOBEC == TRUE)
      ) %>%
      mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
      filter(apobecratio < 0.5) %>%
      select(hotspot) 
    indel_5M <- indel_5M %>%
      mutate(Apobecratio_lesshalf = hotspot %in% hotspot_APOBEC_lesshalf$hotspot)}
  filename <- file.path(folder_path,"indel_5M.txt")
  write.table(indel_5M,file =filename,row.names = F,quote = F,sep = "\t")
}

#4-6 indel,within注释元件取对应基因，others取距离<10M基因
if(1){
  #最终得到数据框indel_10M
  #hotspot 该hotspot对应的id
  #chrom 
  #start(1-based)
  #end(1-based)
  #pval p value
  #length length of hotspot
  #p.bg Mean background mutation probability
  #k number of mutated samples
  #fdr fdr矫正后的q值
  #element within的元件(hotspot within元件)
  #type 
  #target_gene 假定的靶基因,within注释元件即为对应基因，others为距离10M内的某个基因
  #raw_expr_p 突变组相较于非突变组的wilconxin p
  #raw_expr_log2FC
  #cnnormal_expr_p 突变组相较于非突变组的wilconxin p,仅在拷贝数正常的样本中进行比较
  #cnnormal_expr_log2FC
  #raw_expr_allgene_q 在所有基因(>50%样本中表达)中进行fdr 矫正
  #raw_expr_driveronly_q 仅在drivergene(>50%样本中表达)中进行fdr 矫正
  #cnnormal_expr_allgene_q
  #cnnormal_expr_driveronly_q
  #IRratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于一段回文序列里
  #Mapratio_overhalf 布尔值，该hotspot是否超过百分之50突变位于mappable的区域(GEM mappabilty>0.3 且不在encode(DAC) blacklist regions)
  #Apobecratio_lesshalf 布尔值，该hotspot是否小于百分之50突变来自于APOBEC突变过程
  
  #4-6-1 添加overlapping element,type和注释元件的target_gene
  if(1){
    #0-based ，左闭右开,左+1,右不变
    gr_ele <- GRanges(seqnames = all_elements$V1,ranges = IRanges(start = all_elements$V2+1,end = all_elements$V3))
    #1-based，左闭右闭，左闭右闭
    gr_indelhotspots <- GRanges(seqnames = indel_hotspots_merged$chrom,ranges = IRanges(start = indel_hotspots_merged$start,end = indel_hotspots_merged$end))
    #做overlap 前面的是query,后面的是subject,query within subject
    hits_ele_indelhotspots <- findOverlaps(gr_indelhotspots,gr_ele,type="within")
    #由于一个hotspot可对应多个element，需要整合一下(1个hotspot的element弄到一行)
    overlap <- data.frame(element=all_elements$V4[subjectHits(hits_ele_indelhotspots)],hotspot=rownames(indel_hotspots_merged)[queryHits(hits_ele_indelhotspots)])
    overlap <- overlap %>%
      #添加type
      mutate(type = case_when(
        str_detect(element, "gc46_pc.cds") ~ "CDS",
        str_detect(element, "gc46_pc.ss_cds") ~ "SS_CDS",
        str_detect(element, "gc46_pc.5utr") ~ "UTR5",
        str_detect(element, "gc46_pc.3utr") ~ "UTR3",
        str_detect(element, "gc46.promoter") ~ "Promoter_CDS",
        str_detect(element, "lncrna_exons") ~ "LncRNA",
        str_detect(element, "ss_lncrna") ~ "SS_LncRNA",
        str_detect(element, "lncrna.promoter") ~ "Promoter_LncRNA",
        str_detect(element, "smallrna") ~ "SmallRNA",
        str_detect(element, "mirna") ~ "miRNA",
        str_detect(element, "enhancers") ~ "Enhancer",
        TRUE ~ "Others"  # 其他情况赋值为 Others
      )) %>%
      #添加target_gene
      mutate(target_gene = ifelse(!(type %in% c("Others","Enhancer")),sapply(overlap$element,function(x) strsplit(x,"::")[[1]][3]),
                                  ifelse(type == "Enhancer",sapply(overlap$element,function(x) strsplit(x,"::")[[1]][4]),NA))) %>%
      #将每个hotspot对应的多个element整合到一起
      mutate(key1=paste(element,target_gene,type,sep = ":::")) %>%
      select(hotspot,key1) %>%
      group_by(hotspot) %>%
      summarise(key2=paste(key1,collapse = ";"))
    indel_10M <- indel_hotspots_merged %>%
      mutate(hotspot=rownames(indel_hotspots_merged)) %>%
      left_join(overlap,by="hotspot") %>%
      separate_rows(key2,sep = ";") %>%
      separate(key2,into = c("element","target_gene","type"),sep = ":::") %>%
      mutate(type=ifelse(is.na(type),"Others",type))
  }
  #4-6-2 添加others的target_gene
  if(1){
    #0-based 左闭右开，左+1，右不变
    gr_gene <- GRanges(seqnames = all_gene_interval$chr,ranges = IRanges(start = all_gene_interval$start+1,end = all_gene_interval$end))
    #1-based 左闭右闭，左不变，右不变
    gr_indel_10M <- GRanges(seqnames = indel_10M$chrom,ranges = IRanges(start = indel_10M$start,end = indel_10M$end))
    hits_gene_indel10M <- findOverlaps(query = gr_gene,subject = gr_indel_10M,maxgap = 10000000)
    #添加非pro/utr的target_gene
    #每个非pro/utr的hotspot对应多个target_gene,先把它们整理成每个hotspot一行
    hotspot_target <- data.frame(hotspot=indel_10M$hotspot[subjectHits(hits_gene_indel10M)],others_target=all_gene_interval$gene[queryHits(hits_gene_indel10M)])
    hotspot_target <- hotspot_target %>%
      group_by(hotspot) %>%
      summarise(others_targets=paste(others_target,collapse = ";"))
    indel_10M <- indel_10M %>%
      left_join(hotspot_target,by="hotspot") %>%
      mutate(target_gene=ifelse(is.na(target_gene),others_targets,target_gene)) %>%
      separate_rows(target_gene,sep = ";")
    #只保留target_gene在大于50%样本中表达的行
    if(1){
      # 计算每行值为0的列的数量
      zero_count <- apply(expr == 0, 1, sum)
      # 计算总列数
      total_columns <- ncol(expr)
      # 筛选出大于50%列值为0的行的行名
      rownames  <- rownames(expr[zero_count < total_columns / 2, ])
      # 使用 strsplit 按照 "|" 分割
      split_vec <- strsplit(rownames, "\\|")
      # 扁平化为一个向量
      flattened_vec <- unlist(split_vec)
      indel_10M <- indel_10M %>%
        filter(target_gene %in% flattened_vec)
    }
  }
  #4-6-3 添加raw_expr_p,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_log2FC
  #添加raw_expr_p,raw_expr_log2FC
  if(1){
    results <- rep(NA,nrow(indel_10M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(indel_10M), .combine = c) %dopar% {
      library(dplyr)
      p_expr(indel_10M$hotspot[i], indel_10M$target_gene[i], his, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    indel_10M$p_expr <- results
    indel_10M <- indel_10M %>%
      separate(p_expr,into=c("raw_expr_p","raw_expr_log2FC"),sep = ":")
  }
  #添加cnnormal_expr_p,cnnormal_expr_log2FC
  if(1){
    results <- rep(NA,nrow(indel_10M))
    cl <- makeCluster(40)
    # 注册并行后端
    registerDoParallel(cl)
    # 并行执行 for 循环
    results <- foreach(i = 1:nrow(indel_10M), .combine = c) %dopar% {
      library(dplyr)
      p_expr_cnnormal(indel_10M$hotspot[i], indel_10M$target_gene[i], his, expr_rownames, expr_hugosymbol)
    }
    # 停止并行计算
    stopCluster(cl)
    indel_10M$p_expr_cnnormal <- results
    indel_10M <- indel_10M %>%
      separate(p_expr_cnnormal,into=c("cnnormal_expr_p","cnnormal_expr_log2FC"),sep = ":")
  }
  #添加#raw_expr_allgene_q,raw_expr_driveronly_q,cnnormal_expr_allgene_q,cnnormal_expr_driveronly_q
  indel_10M <- indel_10M %>%
    mutate(raw_expr_allgene_q=p.adjust(raw_expr_p,method = "BH")) %>%
    mutate(raw_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,raw_expr_p,NA)) %>%
    mutate(raw_expr_driveronly_q=p.adjust(raw_expr_driveronly_p,method = "BH")) %>%
    mutate(cnnormal_expr_allgene_q=p.adjust(cnnormal_expr_p,method = "BH")) %>%
    mutate(cnnormal_expr_driveronly_p=ifelse(target_gene %in% all_drivergenes_list,cnnormal_expr_p,NA)) %>%
    mutate(cnnormal_expr_driveronly_q=p.adjust(cnnormal_expr_driveronly_p,method = "BH")) %>%
    select(-raw_expr_driveronly_p,-cnnormal_expr_driveronly_p,-others_targets)
  #4-6-4添加IRratio_overhalf,Mapratio_overhalf,Apobecratio_lesshalf
  if(1){
    #IRratio_overhalf
    hotspot_IRratio_overhalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        IR_muts_n = sum(In_IR == TRUE)
      ) %>%
      mutate(IRratio = IR_muts_n / total_muts_n) %>%
      filter(IRratio >= 0.5) %>%
      select(hotspot) 
    indel_10M <- indel_10M %>%
      mutate(IRratio_overhalf = hotspot %in% hotspot_IRratio_overhalf$hotspot)
    #Mapratio_overhalf
    hotspot_Mapratio_overhalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        mappable_muts_n = sum(Mappable == TRUE)
      ) %>%
      mutate(Mapratio = mappable_muts_n / total_muts_n) %>%
      filter(Mapratio >= 0.5) %>%
      select(hotspot) 
    indel_10M <- indel_10M %>%
      mutate(Mapratio_overhalf = hotspot %in% hotspot_Mapratio_overhalf$hotspot)
    #Apobecratio_lesshalf
    hotspot_APOBEC_lesshalf <- his %>%
      group_by(hotspot) %>%
      summarise(
        total_muts_n = n(),
        apobec_muts_n = sum(APOBEC == TRUE)
      ) %>%
      mutate(apobecratio = apobec_muts_n / total_muts_n) %>%
      filter(apobecratio < 0.5) %>%
      select(hotspot) 
    indel_10M <- indel_10M %>%
      mutate(Apobecratio_lesshalf = hotspot %in% hotspot_APOBEC_lesshalf$hotspot)}
  filename <- file.path(folder_path,"indel_10M.txt")
  write.table(indel_10M,file =filename,row.names = F,quote = F,sep = "\t")
}

#6 分析绘图
if(1){
#6-1 snv_tad
filepath <- file.path(folder_path,"snv_tad.txt")
snv_tad <- fread(filepath)
#无 mut-expr association q <0.1者
qq(snv_tad$raw_expr_p)
title("snv_tad_raw_expr_p")
qq(snv_tad$cnnormal_expr_p)
title("snv_tad_cnnormal_expr_p")

#6-2 snv_5M
filepath <- file.path(folder_path,"snv_5M.txt")
snv_5M <- fread(filepath)
#无 mut-expr association q <0.1者
qq(snv_5M$raw_expr_p)
title("snv_5M_raw_expr_p")
qq(snv_5M$cnnormal_expr_p)
title("snv_tad_cnnormal_expr_p")

#6-3 snv_10M
filepath <- file.path(folder_path,"snv_10M.txt")
snv_10M <- fread(filepath)
#无 mut-expr association q <0.1者
qq(snv_10M$raw_expr_p)
title("snv_10M_raw_expr_p")
qq(snv_10M$cnnormal_expr_p)
title("snv_10M_cnnormal_expr_p")

#6-4 indel_tad
filepath <- file.path(folder_path,"indel_tad.txt")
indel_tad <- fread(filepath)
#无 mut-expr association q <0.1者
qq(indel_tad$raw_expr_p)
title("indel_tad_raw_expr_p")
qq(indel_tad$cnnormal_expr_p)
title("indel_tad_raw_expr_p")

#6-5 indel_5M
filepath <- file.path(folder_path,"indel_5M.txt")
indel_5M <- fread(filepath)
#无 mut-expr association q <0.1者
qq(indel_5M$raw_expr_p)
title("indel_5M_raw_expr_p")
qq(indel_5M$cnnormal_expr_p)
title("indel_5M_cnnormal_expr_p")

#6-6 indel_10M
filepath <- file.path(folder_path,"indel_10M.txt")
indel_10M <- fread(filepath)
#无 mut-expr association q <0.1者
qq(indel_10M$raw_expr_p)
title("indel_10M_raw_expr_p")
qq(indel_10M$cnnormal_expr_p)
title("indel_10M_cnnormal_expr_p")

#6-7 snv_hotspots_merged
#6-7-1 qqplot
qq(snv_hotspots_merged$fdr)

#6-8 indel_hotspots_merged
qq(indel_hotspots_merged$fdr)
}

#草稿
if(0){
  pro_and_utr <- all_elements %>%
    filter(str_detect(V4,"promoter|3utr|5utr"))
#去掉NA的行(因为不是每个元件上都有突变)
a <- GRanges(seqnames = c("chr1"),ranges = IRanges(start = c(1),end = c(4)))
b <- GRanges(seqnames = c("chr1","chr1"),ranges = IRanges(start = c(1,2),end = c(3,4)))

hss <- na.omit(hss)
  hits <- data.frame(query=queryHits(hits_tad_snvtad),subject=subjectHits(hits_tad_snvtad))
snv_hotspots_filepath <- file.path(folder_path,"986_jiaoji_result/snv_hotspots.tsv")
#2-1 snv hotspots
snv_hotspots <- fread(snv_hotspots_filepath,data.table = F)
rownames(snv_hotspots) <- snv_hotspots$V1
snv_hotspots$V1 <- NULL
duplicated_rows <- cleaned_df %>%
  filter(duplicated(cbind(hotspot, target_gene)))
cleaned_df <- overlap %>%
  mutate(key=paste(hotspot, target_gene,sep = ":")) %>%
  arrange(
    factor(type, levels = c("Promoter_CDS", "Promoter_LncRNA"), ordered = TRUE),
    # 可以根据其他列进一步排序，默认会按升序排列
    .by_group = TRUE
  ) %>%
  filter(!(duplicated(key) & type %in% c("UTR5", "UTR3"))) %>%
  ungroup() %>%
  select(pro_utr:type)

gr_all_elements <- GRanges(seqnames = all_elements$V4,ranges = IRanges(start = all_elements$V2,end = all_elements$V3))

gr_pr <- GRanges(seqnames = pro_and_utr$V4,ranges = IRanges(start =pro_and_utr$V2,end = pro_and_utr$V3))
#相同靶基因，promoter>utr，将utr去掉
mutate(key=paste(hotspot, target_gene,sep = ":")) %>%
  arrange(
    factor(type, levels = c("Promoter_CDS", "Promoter_LncRNA"), ordered = TRUE),
    .by_group = TRUE
  ) %>%
  filter(!(duplicated(key) & type %in% c("UTR5", "UTR3"))) %>%
  ungroup() %>%
  
}