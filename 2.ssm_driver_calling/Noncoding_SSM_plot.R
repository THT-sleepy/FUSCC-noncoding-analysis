#Noncoding ssm 部分绘图

##1导入包和文件路径和自写函数
#install packages
if(0){
Biocductor_packages <- c("qqman","dplyr","ggplot2",
                         "gt","gtExtras","htmltools",
                         "data.table","ggrepel","clusterProfiler",
                         "org.Hs.eg.db","writexl"
)
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}
for (pkg in c(Biocductor_packages)){
  require(pkg,character.only=T) 
}}
#包
if(1){
  library(qqman)
  library(dplyr)
  library(ggplot2)
  library(gt)
  library(gtExtras)
  library(htmltools)
  library(data.table)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db) #人类基因组注释包
  library(writexl)
  draw_expr_cnnormal <- function(gene,mutated_samples,nc_rownames,nc_hugosymbol,cnvfile){
    #所有的比较限于cnnormal的样本中
    #基因需要在cnv文件中有才行
    if(gene %in% rownames(cnvfile)){
      cn_samples <- names(cnvfile)[cnvfile[gene, ] == 0]
      cn_mutated_samples <- intersect(cn_samples,mutated_samples)
      cn_unmutated_samples <- setdiff(cn_samples,mutated_samples)
      #两个样本集不能为空
      if(!(is_empty(cn_mutated_samples) | is_empty(cn_unmutated_samples))){
        # 获取在表达矩阵中对应的列索引
        mutated_indices <- which(names(normalized_counts) %in% cn_mutated_samples)
        unmutated_indices <- which(names(normalized_counts) %in% cn_unmutated_samples)
        #首先匹配到归一化counts矩阵该gene对应的行
        d_f <- data.frame()
        #如果该基因是ENSG开头,与在rownames里即可
        if(grepl("^ENSG",gene)){
          d_f <- normalized_counts[grepl(gene,nc_rownames),]
        }else{ #如果是Hugo_Symbol与Hugo_Symbol比对，需完全匹配
          d_f <- normalized_counts[nc_hugosymbol==gene,]
        }
        #如果能匹配上
        if(nrow(d_f)==1){
          # 将表达量分组
          d_f["mutation_status",mutated_indices] <- "mutated"
          d_f["mutation_status",unmutated_indices] <- "unmutated"
          d_f1 <- t(as.matrix(d_f))
          d_f1 <- na.omit(d_f1)
          d_f1 <- as.data.frame(d_f1)
          names(d_f1)[1] <- "normalized_counts"
          d_f1$normalized_counts <- as.numeric(d_f1$normalized_counts)
          #去掉表达值为0的
          d_f1 <- d_f1[d_f1$normalized_counts != 0,]
        }else{return(NA)}
      }else{return(NA)}
    }else{return(NA)}
    d_f1 <- as.data.frame(d_f1)
    p <- ggplot(data=d_f1,mapping = aes(x=mutation_status,y=normalized_counts))+
      geom_violin(aes(fill=mutation_status))+
      #设置颜色
      scale_fill_manual(values = c("#A3C6E0", "#F2B2A2"))+
      geom_boxplot(width=0.1)+
      geom_jitter(width=0.2,cex=0.4)+
      #选择主题
      theme_bw()+
      #去掉网格线,纵坐标以及图例
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none",
            axis.ticks.y = element_blank())+
      #纵坐标取对数
      scale_y_log10()+
      #设置坐标轴标题和主标题
      labs(x="",y="",title =paste(gene,"(CN=2)"))
    my_comparisons <- list(c("mutated", "unmutated"))
    p1 <- p + stat_compare_means(comparisons=my_comparisons,
                                 method="wilcox.test",
                                 label="p.format",
    )
    filename=paste0(gene,"expr_boxplot.pdf")
    filename <- file.path(folder_path,filename)
    pdf(filename,width = 8,height = 6,onefile = F)
    print(p1)
    dev.off()}
}
#文件路径
if(1){
  folder_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/SSM图片，第二次组会前/"
  enpm_result_filepath <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/enpm_output.tsv"
  normalized_counts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/counts_1000_log2tpmplus1.txt"
  cnv_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_data_by_genes.txt"
  config_filepath <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/wgs_samples_986.txt"
  mutspot_snv_tad_filepath <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/snv_tad.txt"
  mutspot_indel_tad_filepath <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/indel_tad.txt"
  adw_result_filepath <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/adw_results_jiaoji_2025_0119.txt"
}
#文件预处理
if(1){
  #a 得到config文件
  config <- read.delim(config_filepath,header = F)
  config$V1 <- sapply(config$V1, function(x) {str_remove(x, "(LC).*")})
  
  #b 得到归一化后的counts
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
  nc_rownames <- rownames(normalized_counts)
  nc_hugosymbol <- sapply(nc_rownames ,function(x) strsplit(x,"\\|")[[1]][2])
  #c 拷贝数结果
  cnv <- fread(cnv_filepath,data.table = F) #里面的值是log2(CN)-1
  rownames(cnv) <- cnv$`Gene Symbol`
  cnv[,1:3] <- NULL
}

##2绘图
#2-1 enpm
enpm_result <- read.delim(enpm_result_filepath)
if(0){
  #保存成xls
  a <- enpm_result %>%
    filter(!is.na(raw_expr_p))
  filename <- "/home/data/t190513/1000_noncoding/nc_point_mutations/enpm_output.xlsx"
  write_xlsx(a,path=filename)
} #xlsx在enpm脚本中完成了
#2-1-2 高频突变的情况
#2-1-2-1 高频突变表格
my_theme <- function(data) {
  tab_options(
    data = data,
    heading.title.font.size = 16,
    heading.align = "left",
    table.font.size = 12
  )
}
gt_tab <- enpm_result %>%
  mutate(`Location(hg38)`=Location) %>%
  arrange(desc(mut.n)) %>%
  filter(mut.n >11 & !is.na(raw_expr_p) & Type!="Others") %>%
  filter(Mappable==T & IR == F) %>%
  mutate(`raw_expr_q(c)`=candidates_raw_expr_q,`cn_expr_q(c)`=candidates_cnnormal_expr_q,cn_expr_log2fc=cnnormal_expr_log2fc,cn_expr_p=cnnormal_expr_p) %>%
  dplyr::select(`Location(hg38)`,Type,mut.n,Target_Gene,raw_expr_p,`raw_expr_q(c)`,raw_expr_log2fc,cn_expr_p,`cn_expr_q(c)`,cn_expr_log2fc) %>% 
  gt() %>% 
  fmt_number(columns = raw_expr_p:cn_expr_log2fc, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,`raw_expr_q(c)`, cn_expr_p,`cn_expr_q(c)`)) %>%
  tab_header(title = "Recurrent Point Mutations") %>%
  my_theme()
filepath <- file.path(folder_path,"Recurrent_Point_Mutations.png")
#save_html(html,file =filepath )
gt_tab |> gtsave(filepath,expand = 10)#,vwidth = 16000, vheight = 9000,
#2-1-2-1 高频突变靶基因富集分析
if(1){
  #2-1-2-1-1 得到基因集
  df <- enpm_result %>%
    filter(!is.na(raw_expr_p)) %>%
    filter(Mappable==T & IR == F)
  gene_set <- df$Target_Gene
  #因为里面既有SYMBOL又有ENSG，需要拆开来转换再合并起来
  gene_set_symbol <- gene_set[!str_detect(gene_set,"ENSG")]
  gene_set_ensembl <- gene_set[str_detect(gene_set,"ENSG")]
  gene_set_trans1 <- bitr(geneID = gene_set_symbol,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库
  #这里的ENSEMBLE的几个基因还没办法转换，可能太新了
  if(0){
    gene_set_trans2 <- bitr(geneID = gene_set_ensembl,  #感兴趣的基因集
                            fromType= "ENSEMBL",   #输入ID的类型
                            toType=c("ENTREZID"),   #输出ID的类型，可为多个
                            OrgDb="org.Hs.eg.db")  #物种注释数据库 
  }
  #2-1-2-1-2 GO富集分析
  BP <- enrichGO(gene = gene_set_symbol,  #基因列表(转换的ID)
                 keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,  #物种对应的org包
                 ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                 pvalueCutoff = 0.01,  #p值阈值
                 pAdjustMethod = "fdr",  #多重假设检验校正方式
                 minGSSize = 1,   #注释的最小基因集，默认为10
                 maxGSSize = 500,  #注释的最大基因集，默认为500
                 qvalueCutoff = 0.01,  #q值阈值
                 readable = TRUE)  #基因ID转换为基因名
  df <- BP@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"enpm_GO.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
  #2-1-2-1-3 KEGG富集分析
  KEGG <- enrichKEGG(gene = gene_set_trans1$ENTREZID,   #基因列表(同GO) 
                     organism = "hsa",  #物种
                     keyType = "kegg",  #指定的基因ID类型，默认为kegg
                     minGSSize = 1, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.01,  
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.01)
  #绘制柱状图 主要反应p或q
  df <- KEGG@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"enpm_KEGG.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
}
#2-1-3 按照cnnormal_p进行排序
#这里不存在cnnormal_p为NA而raw_p<0.05的情况，所以不用按raw_p排序画
html <- enpm_result %>%
  mutate(`Location(hg38)`=Location) %>%
  arrange(cnnormal_expr_p) %>%
  filter(cnnormal_expr_p < 0.05,!is.na(raw_expr_p)) %>%
  filter(Mappable==T & IR == F) %>%
  dplyr::select(`Location(hg38)`,Type,mut.n,Target_Gene,raw_expr_p,raw_expr_q,raw_expr_log2fc,cnnormal_expr_p,cnnormal_expr_q,cnnormal_expr_log2fc) %>% 
  gt() %>% 
  fmt_number(columns = raw_expr_p:cnnormal_expr_log2fc, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,raw_expr_q, cnnormal_expr_p,cnnormal_expr_q)) %>%
  tab_header(title = "Point Mutations Associated with Expression Change of Target Gene")
filepath <- file.path(folder_path,"Point_Mutations_Associated_with_Expression_Change_of_Target_Gene.html")
save_html(html,file =filepath )
#2-1-4 qqplot
filepath <- file.path(folder_path,"enpm_p_qqplot.pdf")
pdf(filepath,width = 8,height = 6)
qq(enpm_result$cnnormal_expr_p)
p_value <- enpm_result$cnnormal_expr_p
z <- qnorm(p_value/2) 
lambda  <- round(median(z^2, na.rm = TRUE) / 0.456, 3)
title("enpm_p_qqplot,lambda=1.633")
dev.off()

#2-2 mutspot
#2-2-1 高频突变
#2-2-1-1-1 snv hotspot
mutspot_snv_tad <- read.delim(mutspot_snv_tad_filepath)
#保存成xlsx文件
filename <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/snv_tad.xlsx"
mutspot_snv_tad %>% 
  filter(!is.na(raw_expr_p)) %>%
  dplyr::select(-id,hotspot,everything()) %>%
  write_xlsx(path = filename)
#绘制表格
my_theme <- function(data) {
  tab_options(
    data = data,
    heading.title.font.size = 16,
    heading.align = "left",
    table.font.size = 12
  )
}
gt_tbl <- mutspot_snv_tad %>%
  mutate(chrstart=paste(chrom,start,sep = ":")) %>%
  mutate(Location=paste(chrstart,end,sep = "-")) %>%
  mutate(fdr_q=fdr) %>%
  mutate(mut.n=k) %>%
  mutate(raw_expr_q=raw_expr_allgene_q) %>%
  mutate(`cn_expr_q(c)`=candidates_cnnormal_expr_q) %>% 
  mutate(`raw_expr_q(c)`=candidates_raw_expr_q) %>%
  mutate(cn_expr_p=cnnormal_expr_p) %>%
  mutate(cn_expr_log2FC=cnnormal_expr_log2FC) %>%
  filter(!is.na(raw_expr_p) & type!="Others" & mut.n >12) %>%
  filter(Mapratio_overhalf==T,Apobecratio_lesshalf==T,IRratio_overhalf==F) %>%
  dplyr::select(Location,mut.n,fdr_q,type,target_gene,raw_expr_p,`raw_expr_q(c)`,raw_expr_log2FC,cn_expr_p,`cn_expr_q(c)`,cn_expr_log2FC) %>%
  gt() %>% 
  fmt_number(columns = raw_expr_p:cn_expr_log2FC, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,`raw_expr_q(c)`, cn_expr_p,`cn_expr_q(c)`)) %>%
  tab_header(title = "Hotspot in annotated regions") %>%
  my_theme 
filepath <- file.path(folder_path,"Hotspot_in_annotated_regions.png")
gt_tbl |> gtsave(filepath,expand = 10)
#save_html(html,file =filepath )
#2-2-1-1-2 snv hotspot 富集分析
if(1){
  #2-2-1-1-2-1 得到基因集
  df <- mutspot_snv_tad %>%
    filter(!is.na(raw_expr_p) & type != "Others") %>%
    filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T)
  gene_set <- df$target_gene
  #因为里面既有SYMBOL又有ENSG，需要拆开来转换再合并起来
  gene_set_symbol <- gene_set[!str_detect(gene_set,"ENSG")]
  gene_set_ensembl <- gene_set[str_detect(gene_set,"ENSG")]
  gene_set_trans1 <- bitr(geneID = gene_set_symbol,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库
  gene_set_trans2 <- bitr(geneID = gene_set_ensembl,  #感兴趣的基因集
                          fromType= "ENSEMBL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库 
  gene_set_union <- union(gene_set_trans1$ENTREZID,gene_set_trans2$ENTREZID)
  #2-1-2-1-2 GO富集分析
  BP <- enrichGO(gene = gene_set_union,  #基因列表(转换的ID)
                 keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,  #物种对应的org包
                 ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                 pvalueCutoff = 0.01,  #p值阈值
                 pAdjustMethod = "fdr",  #多重假设检验校正方式
                 minGSSize = 1,   #注释的最小基因集，默认为10
                 maxGSSize = 500,  #注释的最大基因集，默认为500
                 qvalueCutoff = 0.01,  #q值阈值
                 readable = TRUE)  #基因ID转换为基因名
  df <- BP@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"mutspot_snv_GO.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
  #2-1-2-1-3 KEGG富集分析
  KEGG <- enrichKEGG(gene = gene_set_union,   #基因列表(同GO) 
                     organism = "hsa",  #物种
                     keyType = "kegg",  #指定的基因ID类型，默认为kegg
                     minGSSize = 1, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.01,  
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.01)
  #绘制柱状图 主要反应p或q
  df <- KEGG@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"mutspot_snv_KEGG.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
}

#2-2-1-2-1 indel hotspot
mutspot_indel_tad <- read.delim(mutspot_indel_tad_filepath)
#保存成xlsx文件
filename <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/indel_tad.xlsx"
mutspot_indel_tad %>% 
  filter(!is.na(raw_expr_p)) %>%
  dplyr::select(-id,hotspot,everything()) %>%
  write_xlsx(path = filename)
#绘制表格
my_theme <- function(data) {
  tab_options(
    data = data,
    heading.title.font.size = 16,
    heading.align = "left",
    table.font.size = 12
  )
}
gt_tbl <- mutspot_indel_tad %>%
  mutate(chrstart=paste(chrom,start,sep = ":")) %>%
  mutate(Location=paste(chrstart,end,sep = "-")) %>%
  mutate(fdr_q=fdr) %>%
  mutate(mut.n=k) %>%
  mutate(`raw_expr_q(c)`=candidates_raw_expr_q) %>%
  mutate(`cn_expr_q(c)`=candidates_cnnormal_expr_q) %>%
  mutate(cn_expr_log2FC=cnnormal_expr_log2FC) %>%
  mutate(cn_expr_p=cnnormal_expr_p) %>%
  filter(!is.na(raw_expr_p) & type!="Others" & mut.n >5) %>% 
  filter(Mapratio_overhalf==T,Apobecratio_lesshalf==T,IRratio_overhalf==F) %>%
  dplyr::select(Location,mut.n,fdr_q,type,target_gene,raw_expr_p,`raw_expr_q(c)`,raw_expr_log2FC,cn_expr_p,`cn_expr_q(c)`,cn_expr_log2FC) %>%
  gt() %>% 
  fmt_number(columns = raw_expr_p:cn_expr_log2FC, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,`raw_expr_q(c)`, cn_expr_p,`cn_expr_q(c)`)) %>%
  tab_header(title = "Hotspot(indel) in annotated regions") %>%
  my_theme()
filepath <- file.path(folder_path,"Hotspot(indel)_in_annotated_regions.png")
gt_tbl |> gtsave(filepath,expand = 10)
#save_html(html,file =filepath )
#2-2-1-2-2 mutspot indel 富集分析
if(1){
  #2-2-1-2-2-1 得到基因集
  df <- mutspot_indel_tad %>%
    filter(!is.na(raw_expr_p) & type != "Others") %>%
    filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T)
  gene_set <- df$target_gene
  #因为里面既有SYMBOL又有ENSG，需要拆开来转换再合并起来
  gene_set_symbol <- gene_set[!str_detect(gene_set,"ENSG")]
  gene_set_ensembl <- gene_set[str_detect(gene_set,"ENSG")]
  gene_set_trans1 <- bitr(geneID = gene_set_symbol,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库
  if(0){
    gene_set_trans2 <- bitr(geneID = gene_set_ensembl,  #感兴趣的基因集
                            fromType= "ENSEMBL",   #输入ID的类型
                            toType=c("ENTREZID"),   #输出ID的类型，可为多个
                            OrgDb="org.Hs.eg.db")  #物种注释数据库 
    gene_set_union <- union(gene_set_trans1$ENTREZID,gene_set_trans2$ENTREZID)
  } #ensembl的转不了，基因太新了
  #2-2-1-2-2-1 GO富集分析
  BP <- enrichGO(gene = gene_set_trans1$ENTREZID,  #基因列表(转换的ID)
                 keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,  #物种对应的org包
                 ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                 pvalueCutoff = 0.01,  #p值阈值
                 pAdjustMethod = "fdr",  #多重假设检验校正方式
                 minGSSize = 1,   #注释的最小基因集，默认为10
                 maxGSSize = 500,  #注释的最大基因集，默认为500
                 qvalueCutoff = 0.01,  #q值阈值
                 readable = TRUE)  #基因ID转换为基因名
  df <- BP@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"mutspot_indel_GO.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
  #2-2-1-2-2-2 KEGG富集分析
  KEGG <- enrichKEGG(gene = gene_set_trans1$ENTREZID,   #基因列表(同GO) 
                     organism = "hsa",  #物种
                     keyType = "kegg",  #指定的基因ID类型，默认为kegg
                     minGSSize = 1, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.01,  
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.01)
  #绘制柱状图 主要反应p或q
  df <- KEGG@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"mutspot_indel_KEGG.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
}
#2-2-2-1 indel 按照cnnormal_p进行排序
mutspot_indel_tad <- read.delim(mutspot_indel_tad_filepath)
html <- mutspot_indel_tad %>%
  mutate(chrstart=paste(chrom,start,sep = ":")) %>%
  mutate(Location=paste(chrstart,end,sep = "-")) %>%
  mutate(fdr_q=fdr) %>%
  mutate(mut.n=k) %>%
  mutate(raw_expr_q=raw_expr_allgene_q) %>%
  mutate(cnnormal_expr_q=cnnormal_expr_allgene_q) %>% 
  arrange(cnnormal_expr_p) %>%
  filter(!is.na(raw_expr_p) & cnnormal_expr_p < 0.05 & type != "Others") %>%
  filter(Mapratio_overhalf==T,Apobecratio_lesshalf==T,IRratio_overhalf==F) %>%
  dplyr::select(Location,mut.n,fdr_q,type,target_gene,raw_expr_p,raw_expr_q,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_q,cnnormal_expr_log2FC) %>%
  gt() %>% 
  fmt_number(columns = raw_expr_p:cnnormal_expr_log2FC, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,raw_expr_q, cnnormal_expr_p,cnnormal_expr_q)) %>%
  tab_header(title = "Hotspot(indel) Mutations Associated with Expression Change of Target Gene")
filepath <- file.path(folder_path,"Hotspot(indel)_Mutations_Associated_with_Expression_Change_of_Target_Gene.html")
save_html(html,file =filepath )
#基因太少，不用富集分析

#2-2-2-2 snv 按照cnnormal_p进行排序
mutspot_snv_tad <- read.delim(mutspot_snv_tad_filepath)
html <- mutspot_snv_tad %>%
  mutate(chrstart=paste(chrom,start,sep = ":")) %>%
  mutate(Location=paste(chrstart,end,sep = "-")) %>%
  mutate(fdr_q=fdr) %>%
  mutate(mut.n=k) %>%
  mutate(raw_expr_q=raw_expr_allgene_q) %>%
  mutate(cnnormal_expr_q=cnnormal_expr_allgene_q) %>% 
  arrange(cnnormal_expr_p) %>%
  filter(!is.na(raw_expr_p) & cnnormal_expr_p < 0.05 & type != "Others") %>%
  filter(Mapratio_overhalf==T,Apobecratio_lesshalf==T,IRratio_overhalf==F) %>%
  dplyr::select(Location,mut.n,fdr_q,type,target_gene,raw_expr_p,raw_expr_q,raw_expr_log2FC,cnnormal_expr_p,cnnormal_expr_q,cnnormal_expr_log2FC) %>%
  gt() %>% 
  fmt_number(columns = raw_expr_p:cnnormal_expr_log2FC, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,raw_expr_q, cnnormal_expr_p,cnnormal_expr_q)) %>%
  tab_header(title = "Hotspot Mutations Associated with Expression Change of Target Gene")
filepath <- file.path(folder_path,"Hotspot_Mutations_Associated_with_Expression_Change_of_Target_Gene.html")
save_html(html,file =filepath )
#相应的富集分析
if(1){
  # 得到基因集
  df <- mutspot_snv_tad %>%
    filter(!is.na(raw_expr_p) & type != "Others" & cnnormal_expr_p < 0.05) %>%
    filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T)
  gene_set <- df$target_gene
  #因为里面既有SYMBOL又有ENSG，需要拆开来转换再合并起来
  gene_set_symbol <- gene_set[!str_detect(gene_set,"ENSG")]
  gene_set_ensembl <- gene_set[str_detect(gene_set,"ENSG")]
  gene_set_trans1 <- bitr(geneID = gene_set_symbol,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库
  if(0){
    gene_set_trans2 <- bitr(geneID = gene_set_ensembl,  #感兴趣的基因集
                            fromType= "ENSEMBL",   #输入ID的类型
                            toType=c("ENTREZID"),   #输出ID的类型，可为多个
                            OrgDb="org.Hs.eg.db")  #物种注释数据库 
    gene_set_union <- union(gene_set_trans1$ENTREZID,gene_set_trans2$ENTREZID)
  } #ensembl的转不了，基因太新了
  #2-2-1-2-2-1 GO富集分析
  BP <- enrichGO(gene = gene_set_trans1$ENTREZID,  #基因列表(转换的ID)
                 keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,  #物种对应的org包
                 ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                 pvalueCutoff = 0.01,  #p值阈值
                 pAdjustMethod = "fdr",  #多重假设检验校正方式
                 minGSSize = 1,   #注释的最小基因集，默认为10
                 maxGSSize = 500,  #注释的最大基因集，默认为500
                 qvalueCutoff = 0.01,  #q值阈值
                 readable = TRUE)  #基因ID转换为基因名
  df <- BP@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"mutspot_snvexpr_GO.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
  #2-2-1-2-2-2 KEGG富集分析
  KEGG <- enrichKEGG(gene = gene_set_trans1$ENTREZID,   #基因列表(同GO) 
                     organism = "hsa",  #物种
                     keyType = "kegg",  #指定的基因ID类型，默认为kegg
                     minGSSize = 1, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.01,  
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.01)
  #绘制柱状图 主要反应p或q
  df <- KEGG@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"mutspot_snvexpr_KEGG.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
}


#5M和10M在文章里提一嘴就行了，没有新的发现

#2-2-3 pvalue qqplot
#2-2-3-1 mut p
filepath <- file.path(folder_path,"mutspot_hotspot_p_qqplot.pdf")
pdf(filepath,width = 8,height = 6)
qq(unique(mutspot_snv_tad$pval))
p_value <- unique(mutspot_snv_tad$pval)
z <- qnorm(p_value/2) 
lambda  <- round(median(z^2, na.rm = TRUE) / 0.456, 3)
title("mutspot_hotspot_p_qqplot,lambda=76.051")
dev.off()
#2-2-3-2 mut-expr p
filepath <- file.path(folder_path,"mutspot_mutexpr_p_qqplot.pdf")
pdf(filepath,width = 8,height = 6)
qq(unique(mutspot_snv_tad$cnnormal_expr_p))
p_value <- unique(mutspot_snv_tad$cnnormal_expr_p)
z <- qnorm(p_value/2) 
lambda  <- round(median(z^2, na.rm = TRUE) / 0.456, 3)
title("mutspot_mutexpr_p_qqplot,lambda=1.156")
dev.off()

#2-3 activedriverwgs
adw_results <- read.delim(adw_result_filepath)
#保存成xlsx文件
filename <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/1-1 高频非编码点突变分析/adw_result.xlsx"
adw_results %>%
  filter(!is.na(raw_expr_p)) %>%
  write_xlsx(path=filename)

#2-3-1 高频突变元件
#2-3-1-1 promoter
df <- adw_results %>% 
  filter(Type %in% c("Promoter_CDS","Promoter_LncRNA") & Post_filter_q < 0.01 & !is.na(raw_expr_p)) %>%
  filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T)
df_label <- df[df$Post_filter_q <1e-10, ]
filename <- file.path(folder_path,"adw_promoter.pdf")
pdf(filename,width = 8,height = 6)
ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
  geom_point()+
  theme_bw()+ #选择主题
  theme(panel.grid = element_blank())+ #去掉网格线
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.position = c(0.2,0.9)
  )+
  labs(x="")+
  scale_color_manual(values = c("Promoter_CDS" = "blue", "Promoter_LncRNA" = "orange")) +
  geom_label(data = df_label, aes(label = Target_gene), vjust = 1.5, hjust = 1,size=4,show.legend = F)
dev.off()

#2-3-1-2 utr
df <- adw_results %>% 
  filter(Type %in% c("UTR5","UTR3") & Post_filter_q < 0.01  & !is.na(raw_expr_p))%>%
  filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T) 
df_label <- df[df$Post_filter_q <1e-10, ]
filename <- file.path(folder_path,"adw_utr.pdf")
pdf(filename,width = 8,height = 6)
ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
  geom_point()+
  theme_bw()+ #选择主题
  theme(panel.grid = element_blank())+ #去掉网格线
  scale_color_manual(values = c("UTR5" = "blue", "UTR3" = "orange")) +
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.position = c(0.2,0.9)
  )+
  labs(x="")+
  geom_label(data = df_label, aes(label = Target_gene), vjust = 1.5, hjust = 1,size=4,show.legend = F)
dev.off()
#2-3-1-3 SS
df <- adw_results %>% 
  filter(Type %in% c("SS_CDS","SS_LncRNA") & Post_filter_q < 0.01 & !is.na(raw_expr_p)) %>%
  filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T) 
df_label <- df[df$Post_filter_q <1e-10, ]
filename <- file.path(folder_path,"adw_ss.pdf")
pdf(filename,width = 8,height = 6)
ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
  geom_point()+
  theme_bw()+ #选择主题
  theme(panel.grid = element_blank())+ #去掉网格线
  scale_color_manual(values = c("SS_CDS" = "blue", "SS_LncRNA" = "orange")) +
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.position = c(0.2,0.9)
  )+
  labs(x="")+
  geom_label(data = df_label, aes(label = Target_gene), vjust = 1.5, hjust = 1,size=4,show.legend = F)
dev.off()
#2-3-1-4 Enhancer
df <- adw_results %>% 
  filter(Type %in% c("Enhancer") & Post_filter_q < 0.01 & !is.na(raw_expr_p)) %>%
  filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T) 
df_label <- df[df$Post_filter_q <1e-10, ]
filename <- file.path(folder_path,"adw_enhancer.pdf")
pdf(filename,width = 8,height = 6)
ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
  geom_point()+
  theme_bw()+ #选择主题
  theme(panel.grid = element_blank())+ #去掉网格线
  scale_color_manual(values = c("Enhancer" = "blue")) +
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.position = c(0.2,0.9)
  )+
  labs(x="")+
  geom_label_repel(data = df_label, aes(label = Target_gene), vjust = 0.5, hjust = 0.3,size=3,show.legend = F,max.overlaps = 100)
dev.off()
#2-3-1-5 LncRNA
df <- adw_results %>% 
  filter(Type %in% c("LncRNA") & Post_filter_q < 0.01 & !is.na(raw_expr_p)) %>%
  filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T) 
df_label <- df[df$Post_filter_q <1e-8, ]
filename <- file.path(folder_path,"adw_lncrna.pdf")
pdf(filename,width = 8,height = 6)
ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
  geom_point()+
  theme_bw()+ #选择主题
  theme(panel.grid = element_blank())+ #去掉网格线
  scale_color_manual(values = c("LncRNA" = "blue")) +
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.position = c(0.2,0.9)
  )+
  labs(x="")+
  geom_label(data = df_label, aes(label = Target_gene), vjust = 1.5, hjust = 1,size=4,show.legend = F)
dev.off()
#2-3-1-6 miRNA
df <- adw_results %>% 
  filter(Type %in% c("miRNA") & Post_filter_q < 0.01 & !is.na(raw_expr_p)) %>%
  filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T) 
df_label <- df[df$Post_filter_q <1e-3, ]
filename <- file.path(folder_path,"adw_mirna.pdf")
pdf(filename,width = 8,height = 6)
ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
  geom_point()+
  theme_bw()+ #选择主题
  theme(panel.grid = element_blank())+ #去掉网格线
  scale_color_manual(values = c("miRNA" = "blue")) +
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.position = c(0.2,0.9)
  )+
  labs(x="")+
  geom_label(data = df_label, aes(label = Target_gene), vjust = 1.5, hjust = 1,size=4,show.legend = F)
dev.off()

#2-3-1-7 smallRNA 无Post_filter_q小于 0.01者
if(0){
  df <- adw_results %>% 
    filter(Type %in% c("smallRNA") & Post_filter_q < 0.01) %>%
    filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T) 
  df_label <- df[df$Post_filter_q <1e-3, ]
  filename <- file.path(folder_path,"adw_smallrna.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(data = df,aes(x=element_muts_obs,y=-log10(Post_filter_q),colour = Type))+
    geom_point()+
    theme_bw()+ #选择主题
    theme(panel.grid = element_blank())+ #去掉网格线
    scale_color_manual(values = c("smallRNA" = "#FFCC99")) +
    geom_label(data = df_label, aes(label = Target_gene), vjust = 1.5, hjust = 1,size=4)
  dev.off()
}

#2-3-1-8 mut qqplot lambda算出来为0，不知为何
filepath <- file.path(folder_path,"adw_mut_p_qqplot.pdf")
pdf(filepath,width = 8,height = 6)
qq(adw_results$Post_filter_p)
p_value <- adw_results$Post_filter_p
z <- qnorm(p_value/2) 
lambda  <- round(median(z^2, na.rm = TRUE) / 0.456, 3)
title("adw_mut_p_qqplot")
dev.off()

#富集分析
if(1){
  df <- adw_results %>% 
    filter(Post_filter_q < 0.01 &! is.na(raw_expr_p)) %>%
    filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T)
  gene_set <- df$Target_gene
  #因为里面既有SYMBOL又有ENSG，需要拆开来转换再合并起来
  gene_set_symbol <- gene_set[!str_detect(gene_set,"ENSG")]
  gene_set_ensembl <- gene_set[str_detect(gene_set,"ENSG")]
  gene_set_trans1 <- bitr(geneID = gene_set_symbol,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库
  gene_set_trans2 <- bitr(geneID = gene_set_ensembl,  #感兴趣的基因集
                          fromType= "ENSEMBL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库 
  gene_set_union <- union(gene_set_trans1$ENTREZID,gene_set_trans2$ENTREZID)
  #2-2-1-2-2-1 GO富集分析
  BP <- enrichGO(gene = gene_set_union,  #基因列表(转换的ID)
                 keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,  #物种对应的org包
                 ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                 pvalueCutoff = 0.01,  #p值阈值
                 pAdjustMethod = "fdr",  #多重假设检验校正方式
                 minGSSize = 1,   #注释的最小基因集，默认为10
                 maxGSSize = 500,  #注释的最大基因集，默认为500
                 qvalueCutoff = 0.01,  #q值阈值
                 readable = TRUE)  #基因ID转换为基因名
  df <- BP@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"adw_GO.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
  #2-2-1-2-2-2 KEGG富集分析
  KEGG <- enrichKEGG(gene = gene_set_union,   #基因列表(同GO) 
                     organism = "hsa",  #物种
                     keyType = "kegg",  #指定的基因ID类型，默认为kegg
                     minGSSize = 1, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.01,  
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.01)
  #绘制柱状图 主要反应p或q
  df <- KEGG@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"adw_KEGG.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
}


#2-3-2 按cnnormal_expr_p排序
my_theme <- function(data) {
  tab_options(
    data = data,
    heading.title.font.size = 16,
    heading.align = "left",
    table.font.size = 8
  )
}
gt_tbl <- adw_results %>%
  mutate(mut.n=element_muts_obs) %>%
  arrange(candidates_cnnormal_expr_q) %>%
  head() %>%
  filter(Mapratio_overhalf==T,IRratio_overhalf==F,Apobecratio_lesshalf==T) %>%
  mutate(`raw_expr_q(c)`=candidates_raw_expr_q) %>%
  mutate(`cn_expr_q(c)`=candidates_cnnormal_expr_q) %>%
  mutate(cn_expr_p=cnnormal_expr_p) %>%
  mutate(cn_expr_log2FC =cnnormal_expr_log2FC) %>%
  dplyr::select(Type,mut.n,Target_gene,raw_expr_p,`raw_expr_q(c)`,raw_expr_log2FC,cn_expr_p,`cn_expr_q(c)`,cn_expr_log2FC) %>% 
  gt() %>% 
  fmt_number(columns = raw_expr_p:cn_expr_log2FC, n_sigfig= 3) %>% 
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,`raw_expr_q(c)`, cn_expr_p,`cn_expr_q(c)`)) %>%
  tab_header(title = "Highly Mutated Elements Associated with Expression Change of Target Gene") %>%
  my_theme()
filepath <- file.path(folder_path,"Highly_Mutated_Elements_Associated_with_Expression_Change_of_Target_Gene.png")
gt_tbl |> gtsave(filepath,expand = 10)

#save_html(html,file =filepath )
#相应富集分析
if(1){
  df <- adw_results %>% 
    filter(Post_filter_q < 0.01 &! is.na(raw_expr_p) & cnnormal_expr_p < 0.05) %>%
    filter(IRratio_overhalf==F & Mapratio_overhalf==T & Apobecratio_lesshalf==T)
  gene_set <- df$Target_gene
  #因为里面既有SYMBOL又有ENSG，需要拆开来转换再合并起来
  gene_set_symbol <- gene_set[!str_detect(gene_set,"ENSG")]
  gene_set_ensembl <- gene_set[str_detect(gene_set,"ENSG")]
  gene_set_trans1 <- bitr(geneID = gene_set_symbol,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库
  gene_set_trans2 <- bitr(geneID = gene_set_ensembl,  #感兴趣的基因集
                          fromType= "ENSEMBL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Hs.eg.db")  #物种注释数据库 
  gene_set_union <- union(gene_set_trans1$ENTREZID,gene_set_trans2$ENTREZID)
  #2-2-1-2-2-1 GO富集分析
  BP <- enrichGO(gene = gene_set_union,  #基因列表(转换的ID)
                 keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,  #物种对应的org包
                 ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                 pvalueCutoff = 0.01,  #p值阈值
                 pAdjustMethod = "fdr",  #多重假设检验校正方式
                 minGSSize = 1,   #注释的最小基因集，默认为10
                 maxGSSize = 500,  #注释的最大基因集，默认为500
                 qvalueCutoff = 0.01,  #q值阈值
                 readable = TRUE)  #基因ID转换为基因名
  df <- BP@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"adw_expr_GO.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
  #2-2-1-2-2-2 KEGG富集分析
  KEGG <- enrichKEGG(gene = gene_set_union,   #基因列表(同GO) 
                     organism = "hsa",  #物种
                     keyType = "kegg",  #指定的基因ID类型，默认为kegg
                     minGSSize = 1, 
                     maxGSSize = 500,
                     pvalueCutoff = 0.01,  
                     pAdjustMethod = "fdr",
                     qvalueCutoff = 0.01)
  #绘制柱状图 主要反应p或q
  df <- KEGG@result
  df1 <- df %>%
    mutate(`-log10(q)`=-log10(qvalue)) %>%
    arrange(desc(`-log10(q)`)) %>%
    filter(qvalue < 0.1) %>%
    head(20)
  df1$Description <- factor(df1$Description,levels = rev(df1$Description))
  filename <- file.path(folder_path,"adw_expr_KEGG.pdf")
  pdf(filename,width = 8,height = 6)
  ggplot(df1,aes(x=`-log10(q)`,y=Description))+
    geom_bar(stat = "identity",fill="darkorange")+
    theme_bw()+
    labs(y="",x="-log10(FDR)")
  dev.off()
}






#草稿
#箱线图函数
gene <- enpm_result$Target_Gene[42]
mutated_samples <- strsplit(enpm_result$Mutated_Samples[42],",")[[1]]

enpm_result %>%
  mutate(`Location(hg38)`=Location) %>%
  arrange(raw_expr_p) %>%
  filter(raw_expr_p < 0.05 & !is.na(raw_expr_p) & is.na(cnnormal_expr_p)) %>%
  dplyr::select(`Location(hg38)`,Type,mut.n,Target_Gene,raw_expr_p,raw_expr_q,raw_expr_log2fc,cnnormal_expr_p,cnnormal_expr_q,cnnormal_expr_log2fc) %>% 
  gt() %>% 
  fmt_number(columns = raw_expr_p:cnnormal_expr_log2fc, n_sigfig= 3) %>%
  gt_theme_espn() %>% 
  gt_hulk_col_numeric(c(raw_expr_p,raw_expr_q, cnnormal_expr_p,cnnormal_expr_q)) %>%
  tab_header(title = "Point Mutations Associated with Expression Change of Target Gene")


Mutaion_candidates <- c("EEF1A1","ENSG00000287608","ATP6V0C","ENSG00000260272","STMN2",
                        "SNTB2","CXCL6","AP1S1","DAAM1","NOS2","KLRD1","ADCY6-DT",
                        "KPNA3","CBLN3","AP3S2","GGN","TTYH1","BHLHA15","FBXO25","THSD4",
                        "PHF19","TUBG2","CLEC14A","ENSG00000259362",
                        "GDF11","AAMDC","H3C10","NEAT1","MIR4479","PLEK2","WDR4","RNF130")
Positive_genes <- c("EEF1A1","SNTB2","CXCL6","AP1S1","KLRD1","AP3S2","BHLHA15",
                    "FBXO25","THSD4","AAMDC","H3C10","NEAT1","PLEK2","RNF130")
Cox


EPK_VCF <- ems %>%
  filter(str_detect(id,"EEF1A1|PLEK2|KRAS")) %>% 
  mutate(CHROM=chr_mut) %>%
  mutate(POS=start_mut) %>%
  mutate(ID=id) %>%
  mutate(REF=ref_mut) %>%
  mutate(ALT=alt_mut) %>%
  mutate(QUAL=".") %>%
  mutate(FILTER=".") %>%
  mutate(INFO=".") %>%
  select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
write.table(EPK_VCF,file ="/home/data/t190513/1000_noncoding/EPK.vcf",row.names = F,quote = F,sep = "\t")
BHLHA15_VCF <- ems %>%
  filter(start_mut==98211596) %>% 
  mutate(CHROM=chr_mut) %>%
  mutate(POS=start_mut) %>%
  mutate(ID=id) %>%
  mutate(REF=ref_mut) %>%
  mutate(ALT=alt_mut) %>%
  mutate(QUAL=".") %>%
  mutate(FILTER=".") %>%
  mutate(INFO=".") %>%
  head(1) %>%
  select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
write.table(BHLHA15_VCF,file ="/home/data/t190513/1000_noncoding/BHLHA15.vcf",row.names = F,quote = F,sep = "\t")
TH_VCF <- ems %>%
  filter(str_detect(id,"TP53|H3C10")) %>% 
  filter(str_detect(id,"ss_cds")) %>%
  mutate(CHROM=chr_mut) %>%
  mutate(POS=start_mut) %>%
  mutate(ID=id) %>%
  mutate(REF=ref_mut) %>%
  mutate(ALT=alt_mut) %>%
  mutate(QUAL=".") %>%
  mutate(FILTER=".") %>%
  mutate(INFO=".") %>%
  select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
write.table(TH_VCF,file ="/home/data/t190513/1000_noncoding/TH.vcf",row.names = F,quote = F,sep = "\t")
ATKF_BED <- ems %>%
  filter(str_detect(id,"AP1S1|THSD4|KLRD1|FBXO25")) %>% 
  filter(str_detect(id,"3utr")) %>%
  mutate(CHROM=chr_mut) %>%
  mutate(START=start_mut-1) %>%
  mutate(END=end_mut) %>%
  mutate(ID=id) %>%
  mutate(REF=ref_mut) %>%
  mutate(ALT=alt_mut) %>%
  select(CHROM,START,END,ID,REF,ALT)
write.table(ATKF_BED,file ="/home/data/t190513/1000_noncoding/ATKF.bed",row.names = F,quote = F,sep = "\t")
KRAS_BED <- ems %>%
  filter(str_detect(id,"KRAS")) %>% 
  filter(str_detect(id,"promoter")) %>%
  mutate(CHROM=chr_mut) %>%
  mutate(START=start_mut-1) %>%
  mutate(END=end_mut) %>%
  mutate(ID=id) %>%
  mutate(REF=ref_mut) %>%
  mutate(ALT=alt_mut) %>%
  select(CHROM,START,END)
write.table(KRAS_BED,file ="/home/data/t190513/1000_noncoding/KRAS.bed",row.names = F,quote = F,sep = "\t")
TH_BED <- ems %>%
  filter(str_detect(id,"TP53|H3C10")) %>% 
  filter(str_detect(id,"ss_cds")) %>%
  mutate(CHROM=chr_mut) %>%
  mutate(START=start_mut-1) %>%
  mutate(END=end_mut) %>%
  mutate(ID=id) %>%
  mutate(REF=ref_mut) %>%
  mutate(ALT=alt_mut) %>%
  select(CHROM,START,END)
write.table(TH_BED,file ="/home/data/t190513/1000_noncoding/TH.bed3",row.names = F,quote = F,sep = "\t")
if(1){
  category <- rep(c("A", "B"), each = 2)
  p <- c(0.01,0.05,0.1,0.15,0.2,0.99)
  
  qq(p,col = ifelse(category == "A", "blue", "green"))
  text(x = negative_log10_expected_p[1], y = negative_log10_observed_p[1], labels = labels[1], pos = 4, col = "blue")
  negative_log10_observed_p <- -log10(p)
  sorted_p <- sort(p)
  # 计算每个数据点对应的百分位
  percentiles <- sapply(p, function(x) {
    # 获取数据在排序后数组中的位置
    rank <- which(sorted_p == x)
    # 计算百分位数
    percentile <- (rank / length(sorted_p)) 
    return(percentile)
  })
  negative_log10_expected_p <- -log10(qunif(percentiles,min = 0,max = 1))
  labels <- c("A","B","C","D")
}