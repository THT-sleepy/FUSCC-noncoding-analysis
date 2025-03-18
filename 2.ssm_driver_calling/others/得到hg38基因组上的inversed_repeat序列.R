##这个脚本用于得到hg38基因组上的inversed repeat序列(palindromic sequence) 要1h
#目标就是得到一个IR_GRCh38.bed

#按照PCAWG2020 noncoding文章的说法
#最小repeat长度应设为6，中间的间隔序列设为4-8

##1导入相应的包和环境变量
library(detectIR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
##2  用detectIR这个包先弄一下
#测试一下
#test_seq <- "AAAAAATTTTTT"
#test_seq <- "NNNNNNTTTTTT"
#test_result <- detectImperfectIR(seq=test_seq,minStemLen = 6,loopLen = 4,maxMismcNum = 0)
#得到hg38chr1-22,chrX,Y的序列
bsg <- BSgenome.Hsapiens.UCSC.hg38
seq_names <- paste0("chr",1:22)
seq_names[23:24] <- c("chrX","chrY")
#得到hg38chr1-22,X,Y所有的最小重复6bp，环4-8bp的回文序列(IR/Palindromic seq.)
gr_merged <- GRanges()
for(k in seq_names){
chr <- k
seq <- bsg[[chr]]
#最小重复6bp，环4-8bp
for(i in 4:8){
detectIR_result <- detectImperfectIR(seq=seq,minStemLen = 6,loopLen = i,maxMismcNum = 0)
gr <- GRanges(seqnames = chr,ranges = IRanges(start = detectIR_result$startPos,end = detectIR_result$endPos))
gr_merged <- c(gr_merged,gr)
}
}
gr_merged <- reduce(gr_merged)
saveRDS(gr_merged,"/home/data/t190513/1000_noncoding/activedriverwgs/IR_hg38.rds")
sum(width(gr_merged))
#这个结果一共有大约1000万个对象,总长2亿，和PCAWG(hg19)的700万个对象，总长1.3亿差不多
#reduce之后hg19的结果总长大约1.1个亿，估计他们把非1-22，X,Y的也算上了吧，总之我的代码没问题


