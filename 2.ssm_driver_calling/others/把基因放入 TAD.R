#把gene放入TAD
#这个脚本输入TAD和gene的bed文件，返回一个gene放入TAD的文件

#0导入必须的包，设置环境变量
if(1){
library(data.table)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)


TAD_dekker2016_hg19_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/A549_dekker_2016_hg19_TAD.bed" 
TAD_dekker2016_hg38_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/A549_dekker_2016_hg38_TAD.bed" 
all_gene_bed4_hg_19_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/all_gene_hg19.bed4"
all_gene_bed4_hg_38_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/all_gene_hg38.bed4"
result_hg19_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/A549_TAD_withgenes_hg19.bed"
result_hg38_path <- "C:/Users/86152/Desktop/八年制卓越计划/非编码区突变在肿瘤发生发展中的作用课题/源文件/A549_TAD_withgenes_hg38.bed"
}

#1 导入文件
TAD_dekker2016_hg19 <- fread(TAD_dekker2016_hg19_path)
TAD_dekker2016_hg38 <- fread(TAD_dekker2016_hg38_path)
all_gene_bed4_hg_19 <- read.delim(all_gene_bed4_hg_19_path,header = F,sep = " ")
all_gene_bed4_hg_38 <- read.delim(all_gene_bed4_hg_38_path,header = F,sep = " ")

#2 生成结果文件
#2-1 首先把各个bed文件设置成Granges对象(因为bed是左闭右开，所以结尾减去1)
gr_TAD_hg19 <- GRanges(seqnames = TAD_dekker2016_hg19$V1,ranges = IRanges(start = TAD_dekker2016_hg19$V2,end = TAD_dekker2016_hg19$V3-1))
gr_TAD_hg38 <- GRanges(seqnames = TAD_dekker2016_hg38$V1,ranges = IRanges(start = TAD_dekker2016_hg38$V2,end = TAD_dekker2016_hg38$V3-1))
gr_genebed4_hg19 <- GRanges(seqnames = all_gene_bed4_hg_19$V1,ranges = IRanges(start = all_gene_bed4_hg_19$V2,end = all_gene_bed4_hg_19$V3-1))
gr_genebed4_hg38 <- GRanges(seqnames = all_gene_bed4_hg_38$V1,ranges = IRanges(start = all_gene_bed4_hg_38$V2,end = all_gene_bed4_hg_38$V3-1))

#2-2 用FindOverlaps函数把基因放入TAD里
hit_hg19 <- findOverlaps(gr_genebed4_hg19,gr_TAD_hg19,type = "within")
hit_hg38 <- findOverlaps(gr_genebed4_hg38,gr_TAD_hg38,type = "within")
genebed4_hg19_rownumbers <- queryHits(hit_hg19)
TAD_hg19_rownumbers <- subjectHits(hit_hg19)
genebed4_hg38_rownumbers <- queryHits(hit_hg38)
TAD_hg38_rownumbers <- subjectHits(hit_hg38)
#生成 hg19的结果 其实就是给TAD文件加一列genes within 
#设置一些中间向量减少内存使用
result_hg19 <- rep(NA,nrow(TAD_dekker2016_hg19))
genename <- all_gene_bed4_hg_19$V4
for(i in 1:length(hit_hg19)){
  #该hit result_hg19目前的值，也就是genes within目前的值 
  a <- result_hg19[TAD_hg19_rownumbers[i]]
  #该hit 对应genename的值
  b <- genename[genebed4_hg19_rownumbers[i]]
  #如果还没赋过值，直接赋值
  if(is.na(a)){
    result_hg19[TAD_hg19_rownumbers[i]] <- b}
  #否则加在后面
  else{
    result_hg19[TAD_hg19_rownumbers[i]] <- paste(a,b,sep = ";")}
}
result_TAD_hg19_bed <- TAD_dekker2016_hg19[,c(-5)]
result_TAD_hg19_bed$genes_within <- result_hg19
colnames(result_TAD_hg19_bed) <- c("chr","start","end","tad_id","genes_within")
write.table(result_TAD_hg19_bed,result_hg19_path,row.names = F,quote = F)

#生成hg38的结果
#设置一些中间向量减少内存使用
result_hg38 <- rep(NA,nrow(TAD_dekker2016_hg38))
genename <- all_gene_bed4_hg_38$V4
for(i in 1:length(hit_hg38)){
  #该hit result_hg19目前的值，也就是genes within目前的值 
  a <- result_hg38[TAD_hg38_rownumbers[i]]
  #该hit 对应genename的值
  b <- genename[genebed4_hg38_rownumbers[i]]
  #如果还没赋过值，直接赋值
  if(is.na(a)){
    result_hg38[TAD_hg38_rownumbers[i]] <- b}
  #否则加在后面
  else{
    result_hg38[TAD_hg38_rownumbers[i]] <- paste(a,b,sep = ";")}
}
result_TAD_hg38_bed <- TAD_dekker2016_hg38[,c(-5)]
result_TAD_hg38_bed$genes_within <- result_hg38
colnames(result_TAD_hg38_bed) <- c("chr","start","end","tad_id","genes_within")
write.table(result_TAD_hg38_bed,result_hg38_path,row.names = F,quote = F)

#草稿
if(1){
test_TAD <- GRanges(seqnames = c("chr1"),ranges = IRanges(start = c(100,205),end = c(200,260)))
test_gene <- GRanges(seqnames = c("chr1","chr1","chr1","chr1","chr1"),
                     ranges = IRanges(start = c(110,150,200,220,99),end = c(150,250,210,250,99)))
test_gene1 <- GRanges(seqnames = c("chr2","chr2","chr2","chr2","chr2"),
                      ranges = IRanges(start = c(110,150,200,220,99),end = c(150,250,210,250,99)))
hit <- findOverlaps(test_gene1,test_TAD,type = "within")
}

ensembleID <- c("ENSG00000290826","ENSG00000238009","ENSG00000239945")
idTable <- bitr(geneID = ensembleID, 
                fromType = "ENSEMBL", 
                toType = c('ENTREZID','SYMBOL','GENETYPE'),
                OrgDb = org.Hs.eg.db)
