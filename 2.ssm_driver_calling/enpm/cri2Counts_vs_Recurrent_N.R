#===============================================================================
#
#          FILE:criCounts_vs_Recurrent_N.R
# 
#         USAGE: ./criCounts_vs_Recurrent_N.R
# 
#   DESCRIPTION: 这个脚本用于绘制所有突变纵轴为counts，横轴为突变数的图
#       OPTIONS: ---
#  REQUIREMENTS: 
#          BUGS: ---
#         NOTES: 目前使用的输入是986 jiaoji hg38
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 1/17/25
#      REVISION:  ---
#      Reference:
#===============================================================================

#0 导入需要的包和函数和路径
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
library(GenomicRanges)}
#1 导入和预处理文件
#need xiugai
cr_input_filepath <- "/home/data/t190513/1000_noncoding/nc_point_mutations/986_jiaoji.cr.input"
annotations_hg38_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_elements_hg38_noctcf.bed"
folder_path <- "/home/data/t190513/1000_noncoding/nc_point_mutations"
#导入注释文件 by THT
annotations_hg38 <- fread(annotations_hg38_filepath)

#设置成Granges对象
gr_anno_hg38 <- GRanges(seqnames = annotations_hg38$V1,ranges = IRanges(start = annotations_hg38$V2+1,end = annotations_hg38$V3))
#处理empi数据框

#导入vep注释文件
vep_result <- fread(cr_input_filepath,data.table = F)
vep_result$location <- paste(vep_result$V2,vep_result$V3,sep = ":")
vep_result <- vep_result[,-c(2,3)]
colnames(vep_result) <- c("Tumor_Sample_Barcode","ref","alt","location")
#去掉indel
vep_result <- vep_result %>%
  filter(nchar(ref) == 1 & nchar(alt) == 1)
#添加自己的注释Anno_by_tht
if(1){
#作为query的注释文件已经弄成在前面弄成granges格式了
#利用vep_result_germline_removed弄一个subject GRanges文件
seqnames <- sapply(strsplit(as.character(vep_result$location), ":"), `[`, 1)
#vep来自vcf，1-based，要减去1
start <- as.integer(sapply(strsplit(as.character(vep_result$location), ":"), `[`, 2))
end <- start
gr_vep_result <- GRanges(seqnames = seqnames,IRanges(start = start,end = end))
#找到hits
hits <- findOverlaps(gr_anno_hg38,gr_vep_result)
anno_rownumbers <- queryHits(hits)
vep_rownumbers <- subjectHits(hits)
#设置一些中间向量减少内存使用
result <- rep(NA,9207560)
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
vep_result$Anno_by_THT <- result}
#添加一列自己的Main_type
if(1){
#根据pcawg的做法，当一个位置对应多个元件时，
#优先级cds>ss_cds>utr5>utr3>pro_cds>lncrna>ss_lncrna>pro_lncrna>small_rna
#micro_rna>enhancer>ctcfbs
vep_result <- vep_result %>%
    mutate(Main_Type = case_when(
      str_detect(Anno_by_THT, "gc46_pc.cds") ~ "CDS",
      str_detect(Anno_by_THT, "gc46_pc.ss_cds") ~ "SS_CDS",
      str_detect(Anno_by_THT, "gc46_pc.5utr") ~ "UTR5",
      str_detect(Anno_by_THT, "gc46_pc.3utr") ~ "UTR3",
      str_detect(Anno_by_THT, "gc46.promoter") ~ "Promoter_CDS",
      str_detect(Anno_by_THT, "lncrna_exons") ~ "LncRNA",
      str_detect(Anno_by_THT, "ss_lncrna") ~ "SS_LncRNA",
      str_detect(Anno_by_THT, "lncrna.promoter") ~ "Promoter_LncRNA",
      str_detect(Anno_by_THT, "smallrna") ~ "SmallRNA",
      str_detect(Anno_by_THT, "mirna") ~ "miRNA",
      str_detect(Anno_by_THT, "enhancers") ~ "Enhancer",
      str_detect(Anno_by_THT, "ctcfbs") ~ "CTCFbs",
      TRUE ~ "Others"  # 其他情况赋值为 Others
    ))  
}

df_forplot <- vep_result[,c(4,6)]
#绘制 count vs n_sample的图 
df_forplot$Coding_or_non_coding <- ifelse(df_forplot$Main_Type=="CDS","coding","non_coding")
counts <- df_forplot %>%
  count(location,Coding_or_non_coding) %>%
  group_by(n,Coding_or_non_coding) %>%
  summarise(count = n())

#前面代码要跑十几二十分钟，节约时间把处理好的表格保存下来
write.table(counts,"counts",row.names = F,quote = F,sep = "\t")
counts <- read.delim("counts")

#这里注意pacthwork包的版本一定不能是1.3,1.3有bug
p1 <- ggplot(counts, aes(x = n, y = count+0.1,fill=Coding_or_non_coding)) +
  geom_bar(stat =  "identity",color="black",size=0.5) +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  )+
  scale_fill_manual(values = c("coding" = "#4A4A4A","non_coding" = "#D3D3D3")) +
  labs(x = "Recurrent_sample_n", y = "count") +
  theme_minimal()
p2<-p1+
    scale_x_break(c(50,250),#截断位置及范围
                     space = 0.1,#间距大小
                     scales = 0.3)#左右显示比例，大于1上面比例大，小于1下面比例大 
filename <- file.path(folder_path,"986_jiaoji_freq.pdf")
pdf(file = filename,width = 9,height = 5,onefile = F)
print(p2)
dev.off()



#草稿
if(0){
  #合并一些作用不大的列，并且弄成一个ID一行 need xiugai
  if(1){
    df_forplot <- vep_result %>%
      pivot_longer(cols = -c(Uploaded_variation, Location, MAX_AF, Main_Type, BIOTYPE, SYMBOL,ID,Anno_by_THT),
                   names_to = "key", values_to = "value") %>%
      group_by(Uploaded_variation,Location, MAX_AF,Main_Type,ID,Anno_by_THT) %>%
      summarise(Other_Columns = paste(unique(value), collapse = ":"), .groups = 'drop')
    write.table(df_forplot,"df_forplot.txt")}
  #去掉germline snp
  #(有rs记录同时没有COSMIC记录且MAX_AF<0.001，有人发现过且不是在癌症患者中发现的突变)
  #这里的rationale是既然筛除panel of normal时是有就去掉，那这里就相当于弄了一个更大的normal群体，直接去掉而不去考虑频率的因素
  vep_result_germline_removed <- vep_result %>%
    filter(!(str_detect(Existing_variation, "rs") & !str_detect(Existing_variation, "COSV")) & MAX_AF < 0.001)  
#一个样本的一个突变可以有多条注释行，所以需要给每个突变设置一个独特的id,便于后续操作
  vep_result <- vep_result %>%
    mutate(ID = as.numeric(factor(paste(Uploaded_variation,Location,Allele, sep = "_"))))
  if(0){
    #检查下共有多少个独特的突变，发现一共有1305638个，是对的
    #length(unique(vep_result$ID))
    #设置一列新的Type用于绘图 5min
    vep_result <- vep_result %>%
      group_by(Location,Allele) %>%
      mutate(Type = case_when(
        any(EXON != "-" & !grepl("UTR", Consequence) & BIOTYPE == "protein_coding") ~ "protein_coding",
        any(grepl("5_prime_UTR_variant", Consequence)) ~ "5UTR",
        any(grepl("3_prime_UTR_variant", Consequence)) ~ "3UTR",
        any(BIOTYPE == "promoter") ~ "promoter",
        any(BIOTYPE == "enhancer") ~ "enhancer",
        any(grepl("RNA", BIOTYPE)) ~ "nc_rna",
        TRUE ~ "others"
      )) %>%
      ungroup()}
a <- GRanges(seqnames = "chr1",
             ranges = IRanges(start = 1,end = 100))
b <- GRanges(seqnames = c("chr1","chr1"),
             ranges = IRanges(start = c(100,50),end = c(100,70)))

max_semicolon_count <- vep_result_germline_removed %>%
  mutate(num_semicolons = str_count(Anno_by_THT, ";"))%>%
  summarise(max_count = max(num_semicolons[! is.na(num_semicolons)])) %>%
  pull(max_count)

p1 <- ggplot(counts[10,43], aes(x = n, y = count+0.1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  )+
  labs(x = "Recurrent_sample_n", y = "count") +
  theme_minimal()
p2<-p1+scale_x_break(c(90,210),#截断位置及范围
                     space = 0.1,#间距大小
                     scales = 0.3)#左右显示比例，大于1上面比例大，小于1下面比例大 

#绘制 non_coding_ratio vs n_sample的图 


#绘制 n_samples >= 30的 sample vs location图-1
if(1){
  mutation_count <- df_forplot %>%
    group_by(Location,Main_Type,Anno_by_THT) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(Ratio = Count / 820)  # 计算比例
  mutation_count$Anno_by_THT[mutation_count$Location=="chr7:55191822"] <- "EGFR L858R"
  ggplot(mutation_count[mutation_count$Count >= 30,], aes(y = reorder(Location, -Count), x = Ratio,fill = Main_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("CDS" = "#3b6291",
                                 "UTR5" = "#943c39",
                                 "UTR3"="#779043",
                                 "SS_CDS"="#bf7334",
                                 "Promoter_CDS" = "#624c7c",
                                 "LncRNA" = "#388498",
                                 "SS_LncRNA"="#82B0D2",
                                 "Promoter_LncRNA"="#BEB8DC",
                                 "SmallRNA"="#E7DAD2",
                                 "miRNA"="#8ECFC9",
                                 "Enhancer" ="#FFBE7A",
                                 "CTCFbs" ="#FA7F6F",
                                 "Others" = "#B0B0B0")) +
    labs(x ="Ratio", y = "Location") +
    geom_text(aes(label = Anno_by_THT), vjust ="middle",hjust="left",size = 3)+
    scale_x_continuous(limits = c(0, 0.3))+  # 设置 x 轴范围
    theme_bw()
}
#绘制 n_samples >= 25的 sample vs location图-2 把Type=others的去掉
if(1){
  mutation_count <- df_forplot %>%
    group_by(Location,Main_Type,Anno_by_THT) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(Ratio = Count / 820)  # 计算比例
  mutation_count$Anno_by_THT[mutation_count$Location=="chr7:55191822"] <- "EGFR L858R"
  ggplot(mutation_count[mutation_count$Count >= 20 & mutation_count$Main_Type != "Others",], aes(y = reorder(Location, -Count), x = Ratio,fill = Main_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("CDS" = "#3b6291",
                                 "UTR5" = "#943c39",
                                 "UTR3"="#779043",
                                 "SS_CDS"="#bf7334",
                                 "Promoter_CDS" = "#624c7c",
                                 "LncRNA" = "#388498",
                                 "SS_LncRNA"="#82B0D2",
                                 "Promoter_LncRNA"="#BEB8DC",
                                 "SmallRNA"="#E7DAD2",
                                 "miRNA"="#8ECFC9",
                                 "Enhancer" ="#FFBE7A",
                                 "CTCFbs" ="#FA7F6F",
                                 "Others" = "#B0B0B0")) +
    labs(x = "Ratio" , y ="Location") +
    geom_text(aes(label = Anno_by_THT), vjust ="middle",hjust="left",size = 3)+
    scale_x_continuous(limits = c(0, 0.5))+  # 设置 x 轴范围
    theme_bw()
}
#绘制 sample vs location图-3 只画和癌基因相关的(目前没有enhancer,Others和CTCFbs)
if(1){
  mutation_count <- df_forplot %>%
    group_by(Location,Main_Type,Anno_by_THT) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(Ratio = Count / 820)  # 计算比例
  mutation_count$Anno_by_THT[mutation_count$Location=="chr7:55191822"] <- "EGFR L858R"
  #中间向量
  a <- mutation_count$Anno_by_THT
  b <- sapply(a,function(x) any(grepl(paste(all_cancer_gene, collapse="|"), x)))
  write.table(b,"b")
  ggplot(mutation_count[b & mutation_count$Count >15,], aes(y = reorder(Location, -Count), x = Ratio,fill = Main_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("CDS" = "#3b6291",
                                 "UTR5" = "#943c39",
                                 "UTR3"="#779043",
                                 "SS_CDS"="#bf7334",
                                 "Promoter_CDS" = "#624c7c",
                                 "LncRNA" = "#388498",
                                 "SS_LncRNA"="#82B0D2",
                                 "Promoter_LncRNA"="#BEB8DC",
                                 "SmallRNA"="#E7DAD2",
                                 "miRNA"="#8ECFC9",
                                 "Enhancer" ="#FFBE7A",
                                 "CTCFbs" ="#FA7F6F",
                                 "Others" = "#B0B0B0")) +
    labs(x = "Ratio" , y ="Location") +
    geom_text(aes(label = Anno_by_THT), vjust ="middle",hjust="left",size = 3)+
    scale_x_continuous(limits = c(0, 0.5))+  # 设置 x 轴范围
    theme_bw()
}}

#服务器中patchwork包有bug，这段代码可以修复
if(0){
guides_build_mod <- function (guides, theme){
  legend.spacing.y <- calc_element("legend.spacing.y", theme)  # modified by me
  legend.spacing.x <- calc_element("legend.spacing.x", theme)  # modified by me
  legend.box.margin <- calc_element("legend.box.margin", theme) %||% 
    margin()
  widths <- exec(unit.c, !!!lapply(guides, gtable_width))
  heights <- exec(unit.c, !!!lapply(guides, gtable_height))
  just <- valid.just(calc_element("legend.box.just", theme))
  xjust <- just[1]
  yjust <- just[2]
  vert <- identical(calc_element("legend.box", theme), "horizontal")
  guides <- lapply(guides, function(g) {
    editGrob(g, vp = viewport(x = xjust, y = yjust, just = c(xjust, 
                                                             yjust), height = if (vert) 
                                                               heightDetails(g)
                              else 1, width = if (!vert) 
                                widthDetails(g)
                              else 1))
  })
  guide_ind <- seq(by = 2, length.out = length(guides))
  sep_ind <- seq(2, by = 2, length.out = length(guides) - 1)
  if (vert) {
    heights <- max(heights)
    if (length(widths) != 1) {
      w <- unit(rep_len(0, length(widths) * 2 - 1), "mm")
      w[guide_ind] <- widths
      w[sep_ind] <- legend.spacing.x
      widths <- w
    }
  }
  else {
    widths <- max(widths)
    if (length(heights) != 1) {
      h <- unit(rep_len(0, length(heights) * 2 - 1), "mm")
      h[guide_ind] <- heights
      h[sep_ind] <- legend.spacing.y
      heights <- h
    }
  }
  write.table(counts,"counts",row.names = F,quote = F,sep = "\t")
  widths <- unit.c(legend.box.margin[4], widths, legend.box.margin[2])
  heights <- unit.c(legend.box.margin[1], heights, legend.box.margin[3])
  guides <- gtable_add_grob(gtable(widths, heights, name = "guide-box"), 
                            guides, t = 1 + if (!vert) 
                              guide_ind
                            else 1, l = 1 + if (vert) 
                              guide_ind
                            else 1, name = "guides")
  gtable_add_grob(guides, element_render(theme, "legend.box.background"), 
                  t = 1, l = 1, b = -1, r = -1, z = -Inf, clip = "off", 
                  name = "legend.box.background")
}
environment(guides_build_mod) <- asNamespace('patchwork')
assignInNamespace("guides_build", guides_build_mod, ns = "patchwork")}