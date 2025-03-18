#===============================================================================
#
#          FILE:enpmi2Explore_noncoding_point_mutations_output.R
# 
#         USAGE:Rscript enpmi2Explore_noncoding_point_mutations_output.R <folderpath>
#               <snv_filename> <outputxlsx_filename> <outputtsv_filename> <config_filename>
# 
#   DESCRIPTION: 这个脚本用于分析点突变及其与靶基因表达的关系
#       OPTIONS: ---
#  REQUIREMENTS: 
#          BUGS: ---
#         NOTES: 这是第二版；删除了不用的自写函数，改变了log2fc的计算方式(用DESeq2标化的reads算,count+了1，函数中加了log2)
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 17/3/25
#      REVISION: ed2
#      Reference:
#===============================================================================

#输入是vep注释后的文件
#0 导入需要的包和自写函数和路径,设置环境变量
if(1)
{
## loading packages
  library(tidyr)
  library(dplyr)
  library(gridExtra)
  library(data.table)
  library(scales)
  library(stringr)
  library(GenomicRanges)
  library(doParallel)
  args <- commandArgs(trailingOnly = TRUE)
## set filepaths 
  #需要自定义的filepaths
  folder_path <-  args[1]
  snv_filepath <- file.path(folder_path,args[2])                           
  outputxlsx_filepath <- file.path(folder_path,args[3]) 
  outputtsv_filepath <- file.path(folder_path,args[4]) 
  config_filepath <- file.path(folder_path,args[5]) 
  #和数据相关的filepaths
  normalized_counts_filepath <- "/home/data/t190513/1000_noncoding/deseq2_normalized_counts.txt"
  clinical_data_file_path <- "/home/data/t190513/1000_noncoding/activedriverwgs/clin_info_1000_with_pathway_with_gene.txt"
  all_singlemuts_sigs_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/SigProfilerAssignment/986_jiaoji_SigProfilerAssignment_results/all_singlemuts_sigs.txt"
  cnv_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_data_by_genes.txt"
  #通用文件filepaths
  annotations_hg38_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_elements_hg38_noctcf.bed"
  TAD_withgenes_hg38_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/A549_TAD_withgenes_hg38.bed"
  all_gene_bed_hg38_file_path <- "/home/data/t190513/1000_noncoding/nc_point_mutations/all_gene_hg38.bed4"
  GEM_mappability_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/GRCh38_mappability_100mer.gem.bed"
  Encode_blacklist_regions_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/encode_hg38.blacklist.bed"
  IR_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/IR_hg38.rds"
  
##define necessary functions
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
          log2_fc <- log2(median(as.numeric(group_mutated))/median(as.numeric(group_unmutated)))
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
              log2_fc <- log2(median(as.numeric(group_mutated))/median(as.numeric(group_unmutated)))
              result[i] <- paste(p,log2_fc,sep = ":")
            }
          }
          return(paste(result,collapse = ";"))
        }else{return(NA)}
      }else{return(NA)}}else{return(NA)}}

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
 
  #导入注释文件 by THT
  annotations_hg38 <- fread(annotations_hg38_filepath)
  #设置成Granges对象(因为bed是左闭右开，0-based,修改为双闭，1-based，需要start+1)
  gr_anno_hg38 <- GRanges(seqnames = annotations_hg38$V1,ranges = IRanges(start = annotations_hg38$V2+1,end = annotations_hg38$V3))
  
  #得到TAD文件及其GRanges对象，修改同上
  TAD_withgenes_hg38 <- fread(TAD_withgenes_hg38_filepath)
  gr_TAD_withgenes_hg38 <- GRanges(seqnames = TAD_withgenes_hg38$chr,ranges = IRanges(start=TAD_withgenes_hg38$start+1,end=TAD_withgenes_hg38$end))
  
  #得到归一化后的counts
  if(1){
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
  #+1
  normalized_counts <- normalized_counts+1
  }
  
  #拷贝数结果
  cnv <- fread(cnv_filepath,data.table = F) #里面的值是log2(CN)-1
  cnv_1 <- cnv
  rownames(cnv_1) <- cnv_1$`Gene Symbol`
  cnv_1[,1:3] <- NULL
  #只留下config有的样本Tumor的列
  mid <- paste0(config$V1,"LC")
  pattern <- paste(mid,collapse = "|")
  cnv_1 <- cnv_1[,grepl(pattern,colnames(cnv_1))]
  
  #Mappability及IR
  GEM_mappability <- fread(GEM_mappability_filepath) #0-based 左开右闭
  names(GEM_mappability) <- c("chr","start","end")
  encode_blacklist_regions <- fread(Encode_blacklist_regions_filepath) #0-based 左开右闭
  names(encode_blacklist_regions) <- c("chr","start","end")
  IR <- readRDS(IR_filepath) #1-based 双闭
  IR <- GenomicRanges::reduce(IR)
}

#===============================================================================

#2 得到输出数据框
#2-1 处理snv数据框，得到表格1 enpm_tad_output
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
#candidates_raw_expr_q 只计算candidates的fdr q值，无拷贝数矫正
#cnnormal_expr_q 突变组与非突变组表达值fdr q值，有拷贝数矫正
#candidates_cnnormal_expr_q 只计算candidates的fdr q值，有拷贝数矫正
#因为这里没有突变和apobec muts重合，所以没有APOBEC列
if(1){
  #2-1-1 导入vep注释文件，包含location
  vep_result <- fread(snv_filepath)
  vep_result$location <- paste(vep_result$V2,vep_result$V3,sep = ":")
  vep_result <- vep_result[,-c(2,3)]
  colnames(vep_result) <- c("Tumor_Sample_Barcode","ref","alt","Location")
  #去掉indel，如果有的话
  vep_result <- vep_result %>%
    filter(nchar(ref) == 1 & nchar(alt) == 1)
  sample_size <- length(unique(vep_result$Tumor_Sample_Barcode))
  #2-1-2 添加Overlap_element
  if(1){
    #作为query的注释文件已经弄成在前面弄成granges格式了
    #利用vep_result弄一个subject GRanges文件
    if(1){
      seqnames <- sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 1)
      #vep来自vcf，1-based，不用改
      start <- as.integer(sapply(strsplit(as.character(vep_result$Location), ":"), `[`, 2))
      end <- start
      gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))}
    #找到hits
    hits <- findOverlaps(gr_anno_hg38,gr_vep_result)
    anno_rownumbers <- queryHits(hits)
    vep_rownumbers <- subjectHits(hits)
    #设置一些中间向量减少内存使用
    result <- rep(NA,nrow(vep_result))
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
  #之所以要分开添加是因为之前id弄的不规范，enhancer的基因在第四个
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
    vep_result$Target_Gene[vep_result$Type == "Others"] <- Others$Target_Gene
    #Others可以有多个Target_gene,拆开成一个一行
    vep_result <- vep_result %>%
      separate_rows(Target_Gene, sep = ";")}
    #去掉mut f小于0.01的行
    n <- ceiling(sample_size * 0.01)
    vep_result <- vep_result[vep_result$mut.n >= n,]  ##args3
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
  #2-1-7-2得到gr_vep_result
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
  vep_result[unique(subjectHits(hits_IR_vep)),"IR"] <- TRUE
  
  #2-1-9 添加靶基因expr的p和log2fc值(without cnv correction) #1h for 1000samples
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
  vep_result$raw_expr_p[vep_result$raw_expr_p %in% c("NaN","NA")] <- NA
 
  
  #2-1-10 添加靶基因expr的p_cnnormal和log2fc值(with cnv correction) #1h for 1000samples
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
  vep_result1$cnnormal_expr_p[vep_result$raw_expr_p %in% c("NaN","NA")] <- NA
  #2-1-11 由于有一些基因在超过50%的基因中表达都为1(实际上是0，加了1)，需要把这些基因算出来的expr_p去掉
  #得到这些基因在表达矩阵中对应的行名
  # 计算每行值为0的列的数量
  zero_count <- apply(normalized_counts == 1, 1, sum)
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
  vep_result1$raw_expr_p <- as.numeric(vep_result1$raw_expr_p)
  vep_result1$raw_expr_q <- p.adjust(vep_result1$raw_expr_p,method = "BH")
  #2-1-13 计算cnnormal expr q值
  vep_result1$cnnormal_expr_p <- as.numeric(vep_result1$cnnormal_expr_p)
  vep_result1$cnnormal_expr_q <- p.adjust(vep_result1$cnnormal_expr_p,method = "BH")
  #2-14 输出结果
  output <- vep_result1 %>%
    filter(!is.na(raw_expr_p))
  output$candidates_raw_expr_q <- NA
  output$candidates_raw_expr_q[output$Mappable==T & output$IR==F] <- p.adjust(output$raw_expr_p[output$Mappable==T & output$IR==F],method = "BH")
  output$candidates_cnnormal_expr_q <- NA
  output$candidates_cnnormal_expr_q[output$Mappable==T & output$IR==F] <- p.adjust(output$cnnormal_expr_p[output$Mappable==T & output$IR==F],method = "BH")
  writexl::write_xlsx(output,path=outputxlsx_filepath)
  write.table(output,outputtsv_filepath,row.names = F,quote = F,sep = "\t")
}
