#===============================================================================
#          FILE:tp53_splice_analysis.R
#         USAGE:
#   DESCRIPTION: 这个脚本用于分析fuscc 1000 cohort里的tp53 splice variants
#  REQUIREMENTS: 编码区突变maf文件;非编码区突变maf文件
#          BUGS: ---
#         NOTES: 
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 12/2/25
#      REVISION:  ---
#      Reference:
#===============================================================================

## loading packages
library(dplyr)
library(pheatmap)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(DESeq2)
library(GSVA)
library('GSEABase')
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(survminer)
library(survival)

## set input filepaths
maf_filepath <- "/home/data/t190513/1000_noncoding/mut_fus_combined.txt"
noncoding_maf_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/ems.tsv"
wgs_986_config_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/wgs_samples_986.txt"
counts_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/counts_1000_log2tpmplus1.txt"
cnv_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/all_data_by_genes.txt"
raw_counts_filepath <- "/home/data/t190513/1000_noncoding/fuscc_lc_counts_2010.txt"
hallmark_tp53_pathway_filepath <- "/home/data/t190513/1000_noncoding/HALLMARK_P53_PATHWAY.v2024.1.Hs.gmt"
clin_data_filepath <- "/home/data/t190513/1000_noncoding/activedriverwgs/clin_info_1000_with_pathway_with_gene.txt"


## set output filepaths
waterfallplot_filepath <- "/home/data/t190513/1000_noncoding/tp53_co_mutation.pdf"
tp53_mutationtype_barplot_filepath <- "/home/data/t190513/1000_noncoding/tp53_mutationtype_barplot.pdf"
tp53_expr_muttype_boxplot_filepath <- "/home/data/t190513/1000_noncoding/tp53_expr_muttype_boxplot.pdf"
vacalo_tp53coding_vs_wild_filepath <- "/home/data/t190513/1000_noncoding/vacalo_tp53coding_vs_wild.pdf"
deseq2_tp53coding_vs_wild.txt_filepath <- "/home/data/t190513/1000_noncoding/1000_deseq2_tp53coding_vs_wild.txt"
deseq2_tp53splice_vs_wild.txt_filepath <- "/home/data/t190513/1000_noncoding/1000_deseq2_tp53splice_vs_wild.txt"
deseq2_tp53splice_vs_coding.txt_filepath <- "/home/data/t190513/1000_noncoding/1000_deseq2_tp53splice_vs_coding.txt"
GSEA_tp53coding_vs_wild_filepath <- "/home/data/t190513/1000_noncoding/GSEA_tp53coding_vs_wild.pdf"
GSEA_tp53splice_vs_wild_filepath <- "/home/data/t190513/1000_noncoding/GSEA_tp53splice_vs_wild.pdf"
GSEA_tp53splice_vs_coding_filepath <- "/home/data/t190513/1000_noncoding/GSEA_tp53splice_vs_coding.pdf"
ssGSEA_tp53_coding_vs_splice_vs_wild_filepath <- "/home/data/t190513/1000_noncoding/ssGSEA_tp53_boxplot.pdf"
tp53_mutationtype_os_km_filepath <- "/home/data/t190513/1000_noncoding/tp53_mutationtype_os_km.pdf"
tp53_mutationtype_rfs_km_filepath <- "/home/data/t190513/1000_noncoding/tp53_mutationtype_rfs_km.pdf"

##loading required files
coding_maf <- read.delim(maf_filepath)
noncoding_maf <- read.delim(noncoding_maf_filepath)

config <- read.delim(wgs_986_config_filepath,header = F)
config$V1 <- sapply(config$V1, function(x) {str_remove(x, "(LC).*")})

expr <- fread(counts_filepath,data.table = F)
rownames(expr) <- expr$V1 #fread不支持行名处理一下
expr$V1 <- NULL #fread不支持行名处理一下
mid <- paste0(config$V1,"_T") #只留下config有的样本Tumor的列
pattern <- paste(mid,collapse = "|") #只留下config有的样本Tumor的列
expr <- expr[,grepl(pattern,colnames(expr))] #只留下config有的样本Tumor的列
colnames(expr) <- sapply(colnames(expr),function(x) paste0(strsplit(x,"_")[[1]][3],strsplit(x,"_")[[1]][2])) #改一下列名

cnv <- fread(cnv_filepath,data.table = F) #里面的值是log2(CN)-1
cnv_1 <- cnv
rownames(cnv_1) <- cnv_1$`Gene Symbol`
cnv_1[,1:3] <- NULL

raw_counts <- fread(raw_counts_filepath,data.table = F)
rownames(raw_counts) <- raw_counts$gene_id
raw_counts$gene_id <- NULL

clin <- read.delim(clin_data_filepath)
clin$sample <- clin$Sample_ID
##get tp53_df
#===============================================================================
#hugo_symbol sample coding_or_splice type
coding_tp53_df <- coding_maf %>%
  filter(Hugo_Symbol == "TP53") %>%
  mutate(hugo_symbol=Hugo_Symbol,sample=Tumor_Sample_Barcode) %>%
  mutate(coding_or_splice="coding",type=Variant_Classification) %>%
  select(hugo_symbol,sample,coding_or_splice,type)
splice_tp53_df <- noncoding_maf %>%
  filter(str_detect(id,"TP53") & str_detect(id,"ss_cds")) %>%
  mutate(hugo_symbol="TP53",sample=sample,coding_or_splice="splice",type="Splice_Mutation") %>%
  select(hugo_symbol,sample,coding_or_splice,type)
tp53_df <- rbind(coding_tp53_df,splice_tp53_df)

##co-mutation analysis(splice vs coding)
#===============================================================================
#plot
pdf(file = waterfallplot_filepath,width = 8,height = 6)
mat <- tp53_df %>%
  group_by(hugo_symbol,sample,coding_or_splice) %>%
  dplyr::count() %>%
  mutate(n=1) %>%
  pivot_wider(names_from = sample,values_from = n,values_fill = list(n = 0))
rownames <- paste(mat$hugo_symbol,mat$coding_or_splice)
mat <- mat[,3:ncol(mat)]
#将突变弄成瀑布样式
new_row <- colSums(mat[,]) * 10
for (i in 1:(ncol(mat))) {
  a <- 0
  for (k in 1:(nrow(mat))){
    if (mat[k,i] != 0){
      a <-  a + 10 - k
    }
  }
  new_row[i] <- new_row[i] + a
}
sorted_columns <- names(mat)[order(-new_row[])]
mat_sorted <- mat[,sorted_columns]
rownames(mat_sorted) <- rownames

pheatmap(as.matrix(mat_sorted),
         cluster_rows = F,
         cluster_cols = F,
         labels_col = rep("",269),
         legend = F,
         cellwidth = 1.2,
         cellheight = 18,
         color = c("#FFCC99","#E5CCFF"))
dev.off()

##type barplot
#===============================================================================
pdf(file = tp53_mutationtype_barplot_filepath,width = 9,height = 5)
df <- tp53_df %>%
    mutate(type = ifelse(is.na(type),"Wild-type",
                       ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"FS",
                              ifelse(type=="Missense_Mutation","MIS",
                                     ifelse(type=="Nonsense_Mutation","NS",
                                            ifelse(type=="Splice_Mutation","SP",
                                                   ifelse(type=="In_Frame_Del","IFD",type))))))) %>%
    select(-coding_or_splice) %>%
    group_by(hugo_symbol,sample) %>%
    summarise(type = paste(type, collapse = ";\n"),
              .groups = "drop")
df$type <- factor(df$type,levels = c("MIS","FS","NS","SP",
                                     "FS;\nSP","MIS;\nMIS","MIS;\nSP",
                                     "FS;\nMIS","IFD","MIS;\nFS",
                                     "MIS;\nFS;\nMIS","MIS;\nIFD","NS;\nMIS",
                                     "NS;\nNS"))
ggplot(df,aes(x=type))+
    geom_bar(fill="#FFCC99",color="black",size=0.2)+
    labs(x="",y="")+
    theme_bw()+
    theme(
    axis.text.x = element_text(size = 15, face = "bold"),  # 设置x轴刻度文字大小和样式
    axis.text.y = element_text(size = 15, face = "bold")   # 设置y轴刻度文字大小和样式
  )
  
dev.off()
##expression box plot(missense frameshift nonsense splice wild-type)
#===============================================================================
#sample type expr(log2(tpm+1)) cn
pdf(file=tp53_expr_muttype_boxplot_filepath,width = 9,height = 6)
tp53_expr <- as.data.frame(t(expr["ENSG00000141510|TP53",]))
tp53_expr$sample <- rownames(tp53_expr)
tp53_cn <- as.data.frame(t(cnv_1["TP53",]))
tp53_cn$sample <- rownames(tp53_cn)
singlehit_samples <- tp53_df %>%
  mutate(type = ifelse(is.na(type),"Wild-type",
                       ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"FS",
                              ifelse(type=="Missense_Mutation","MIS",
                                     ifelse(type=="Nonsense_Mutation","NS",
                                            ifelse(type=="Splice_Mutation","SP",
                                                   ifelse(type=="In_Frame_Del","IFD",type))))))) %>%
  select(-coding_or_splice) %>%
  group_by(hugo_symbol,sample) %>%
  summarise(type = paste(type, collapse = ";"),
            .groups = "drop") %>%
  filter(type %in% c("MIS","FS","NS","SP","IFD")) %>%
  select(sample)
singlehit_samples <- singlehit_samples$sample
df <- config %>%
  mutate(sample = paste0(V1,"LC")) %>%
  left_join(tp53_df,by=c("sample")) %>%
  mutate(type = ifelse(is.na(type),"Wild-type",
      ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"Frame_Shift",
             ifelse(type=="Missense_Mutation","Missense",
                    ifelse(type=="Nonsense_Mutation","Nonsense",
                           ifelse(type=="Splice_Mutation","Splice",type)))))) %>%
  left_join(tp53_expr,by=c("sample")) %>%
  mutate(expr = `ENSG00000141510|TP53`) %>% #log2(TPM+1)
  left_join(tp53_cn,by=c("sample")) %>%
  mutate(cn = TP53) %>% #log2(CN)-1 
  filter(sample %in% singlehit_samples | type=="Wild-type") %>%                            #去除有多个TP53突变的患者
  select(sample,type,expr,cn)
p <- ggplot(data=df,mapping = aes(x=type,y=expr))+
  #geom_violin(aes(fill=type))+
  geom_boxplot(width=0.4,aes(fill=type))+
  #设置颜色
  scale_fill_manual(values = c("#99CCFF", "#FFCC99", "#E5CCFF", "#CCFF99", "#F2B2A2","darkred"))+
  #geom_jitter(width=0.2,cex=0.4)+
  #选择主题
  theme_bw()+
  #去掉网格线,纵坐标以及图例
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.y = element_blank())+
  #设置坐标轴标题和主标题
  labs(x="",y="log2(TPM+1)")+
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  # 设置x轴刻度文字大小和样式
    axis.text.y = element_text(size = 15, face = "bold")   # 设置y轴刻度文字大小和样式
  )
my_comparisons <- list(c("Frame_Shift", "Wild-type"),c("In_Frame_Del", "Wild-type"), 
                       c("Missense", "Wild-type"),c("Nonsense", "Wild-type"),
                       c("Splice","Wild-type")
)
p1 <- p + stat_compare_means(comparisons=my_comparisons,
                             method="wilcox.test",
                             label="p.signif"
)
p1
dev.off()

##Differential expressed gene analysis
#===============================================================================
mid <- paste0(config$V1,"_T") #只留下config有的样本Tumor的列
pattern <- paste(mid,collapse = "|") #只留下config有的样本Tumor的列
t_raw_counts <- raw_counts[,grepl(pattern,colnames(raw_counts))] #只留下config有的样本Tumor的列
colnames(t_raw_counts) <- sapply(colnames(t_raw_counts),function(x) paste0(strsplit(x,"_")[[1]][3],strsplit(x,"_")[[1]][2])) #改一下列名
#去掉既有coding又有splicing的样本
t_raw_counts <- t_raw_counts %>%
  select(-c("1743LC","2173LC","4544LC","4909LC",
           "5373LC","5593LC","6479LC"))
count_matrix <- as.matrix(t_raw_counts)
coldata <- data.frame(sample_type=rep("Tumor",979))
coldata <- coldata %>%
  mutate(sample=colnames(count_matrix)) %>%
  left_join(tp53_df,by=c("sample")) %>% 
  mutate(tp53_type=ifelse(is.na(coding_or_splice),"wild",coding_or_splice)) %>%
  select(sample,sample_type,tp53_type)
coldata <- unique(coldata)
rownames(coldata) <- colnames(count_matrix)

coldata$sample <- NULL
coldata$sample_type <- factor(coldata$sample_type)
coldata$tp53_type <- factor(coldata$tp53_type)

all(colnames(count_matrix) %in% rownames(coldata))#注释文件顺序必须与矩阵的列顺序一致
all(colnames(count_matrix) == rownames(coldata))#注释文件顺序必须与矩阵的列顺序一致
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design= ~ tp53_type)
dds <- DESeq(dds,parallel = T)

#coding vs wild
res <- results(dds,contrast = c('tp53_type', "coding","wild"))
write.table(res,deseq2_tp53coding_vs_wild.txt_filepath,sep="\t", quote=F, col.names=NA) #保留中间结果

df <- read.delim(deseq2_tp53coding_vs_wild.txt_filepath,row.names = 1)
df <- na.omit(df)

cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =0.3 #设置log2FC的阈值
df$Sig = ifelse(df$padj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                  abs(df$log2FoldChange) >= cut_off_log2FC,  #abs绝对值
                ifelse(df$log2FoldChange > cut_off_log2FC ,'Up','Down'),'no')
df$log10FDR <- -log(df$padj, base = 10)
df$log10FDR[is.infinite(df$log10FDR)] <- runif(sum(is.infinite(df$log10FDR)), min = 50, max = 140)
df$symbol <- sapply(rownames(df), function(x) strsplit(x, "\\|")[[1]][2])

#coding vs wild-vacalo plot
pdf(file=vacalo_tp53coding_vs_wild_filepath,width = 8,height = 6) 
p <- ggplot(df, aes(x=log2FoldChange, y=log10FDR, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-10, 10)) +  ylim(0,150)+#调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("TP53 coding vs wild") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()
  )
up_data <- filter(df, Sig == 'Up') %>%  # 从DEG_limma_voom中筛选出上调的基因
  distinct(symbol, .keep_all = TRUE) %>%               # 去除重复的基因，保留第一个出现的行
  top_n(10, log10FDR)                           # 选择-P.Value值最大的前10个基因
down_data <- filter(df, Sig == 'Down') %>%  # 从DEG_limma_voom中筛选出下调的基因
  distinct(symbol, .keep_all = TRUE) %>%                   # 去除重复的基因，保留第一个出现的行
  top_n(10, log10FDR)                               # 选择-P.Value值最大的前10个基因
label <- rbind(up_data,down_data)
p1 <- p +  # 基于普通火山图p
    geom_label_repel(data = label,
                     aes(x=log2FoldChange,y=log10FDR, label = symbol),
                     size = 3, fill="#CCFFFF",
                     box.padding = 0.5, #字到点的距离 
                     min.segment.length = 0.5, #短线段可以省略,
                     max.overlaps=100,
                     segment.color = "black", #segment.colour = NA, 不显示线段
                     show.legend = F)   # 添加上调基因的标签
p1
dev.off()

#coding vs wild-GSEA enrich plot
tp53_geneSets <- read.gmt(hallmark_tp53_pathway_filepath)
df <- read.delim(deseq2_tp53coding_vs_wild.txt_filepath,row.names = 1)
df <- na.omit(df)
df_sort <- df %>%
  mutate(symbol = sapply(rownames(df), function(x) strsplit(x, "\\|")[[1]][2])) %>%
  arrange(desc(log2FoldChange))
gene_list <- df_sort$log2FoldChange
names(gene_list) <- df_sort$symbol
gene_list <- gene_list[gene_list!=0]
gene_list <- gene_list[!duplicated(names(gene_list))]
res <- GSEA(
  gene_list,
  TERM2GENE = tp53_geneSets,
  pvalueCutoff = 0.25
)
pdf(file = GSEA_tp53coding_vs_wild_filepath,width = 8,height = 6)
gseaplot2(res,
          geneSetID = 1,
          title = "Coding vs Wild HALLMARK_TP53_PATHWAY",
          pvalue_table = T,
            )
dev.off()

#coding vs splice vs wild ssGSEA box plot
tp53_geneSets <- read.gmt(hallmark_tp53_pathway_filepath)
df <- expr
df$id <- sapply(rownames(df),function(x) strsplit(x,split = "\\|")[[1]][2])
df <- df[!duplicated(df$id),]
rownames(df) <- df$id
df$id <- NULL
dat <- as.matrix(df)
gsvaPar <- gsvaParam(exprData = dat, 
                     geneSets = tp53_geneSets,
                     kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose = FALSE)
tp53_gsva <- as.data.frame(t(gsva.es))
names(tp53_gsva) <- c("tp53_pathway_enrichment_score")
tp53_gsva$sample <- rownames(tp53_gsva)

singlehit_samples <- tp53_df %>%
  mutate(type = ifelse(is.na(type),"Wild-type",
                       ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"FS",
                              ifelse(type=="Missense_Mutation","MIS",
                                     ifelse(type=="Nonsense_Mutation","NS",
                                            ifelse(type=="Splice_Mutation","SP",
                                                   ifelse(type=="In_Frame_Del","IFD",type))))))) %>%
  select(-coding_or_splice) %>%
  group_by(hugo_symbol,sample) %>%
  summarise(type = paste(type, collapse = ";"),
            .groups = "drop") %>%
  filter(type %in% c("MIS","FS","NS","SP","IFD")) %>%
  select(sample)
singlehit_samples <- singlehit_samples$sample
df <- config %>%
  mutate(sample = paste0(V1,"LC")) %>%
  left_join(tp53_df,by=c("sample")) %>%
  mutate(type = ifelse(is.na(type),"Wild-type",
                       ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"Frame_Shift",
                              ifelse(type=="Missense_Mutation","Missense",
                                     ifelse(type=="Nonsense_Mutation","Nonsense",
                                            ifelse(type=="Splice_Mutation","Splice",type)))))) %>%
  left_join(tp53_gsva,by=c("sample")) %>%
  filter(sample %in% singlehit_samples | type=="Wild-type") %>%#去除有多个TP53突变的患者
  select(sample,type,tp53_pathway_enrichment_score)
p <- ggplot(data=df,mapping = aes(x=type,y=tp53_pathway_enrichment_score))+
  #geom_violin(aes(fill=type))+
  geom_boxplot(width=0.4,aes(fill=type))+
  #设置颜色
  scale_fill_manual(values = c("#99CCFF", "#FFCC99", "#E5CCFF", "#CCFF99", "#F2B2A2","darkred"))+
  #geom_jitter(width=0.2,cex=0.4)+
  #选择主题
  theme_bw()+
  #去掉网格线,纵坐标以及图例
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.y = element_blank())+
  #设置坐标轴标题和主标题
  labs(x="",y="TP53_pathway_enrichment_score")+
  theme(
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),  # 设置x轴刻度文字大小和样式
    axis.text.y = element_text(size = 15, face = "bold")   # 设置y轴刻度文字大小和样式
  )
  
my_comparisons <- list(c("Frame_Shift", "Wild-type"),c("In_Frame_Del", "Wild-type"), 
                       c("Missense", "Wild-type"),c("Nonsense", "Wild-type"),
                       c("Splice","Wild-type"),c("Splice","Frame_Shift"),
                       c("Splice","In_Frame_Del"),
                       c("Splice","Missense"),c("Splice","Nonsense")
                       
)
p1 <- p + stat_compare_means(comparisons=my_comparisons,
                             method="wilcox.test",
                             label="p.signif",
)
pdf(file = ssGSEA_tp53_coding_vs_splice_vs_wild_filepath,width = 9,height = 6)
p1
dev.off()

#splice vs wild
res <- results(dds,contrast = c('tp53_type', "splice","wild"))
write.table(res,deseq2_tp53splice_vs_wild.txt_filepath,sep="\t", quote=F, col.names=NA) #保留中间结果

#splice vs wild-GSEA plot
tp53_geneSets <- read.gmt(hallmark_tp53_pathway_filepath)
df <- read.delim(deseq2_tp53splice_vs_wild.txt_filepath,row.names = 1)
df <- na.omit(df)
df_sort <- df %>%
  mutate(symbol = sapply(rownames(df), function(x) strsplit(x, "\\|")[[1]][2])) %>%
  arrange(desc(log2FoldChange))
gene_list <- df_sort$log2FoldChange
names(gene_list) <- df_sort$symbol
gene_list <- gene_list[gene_list!=0]
gene_list <- gene_list[!duplicated(names(gene_list))]
res <- GSEA(
  gene_list,
  TERM2GENE = tp53_geneSets,
  pvalueCutoff = 0.8
)
pdf(file = GSEA_tp53splice_vs_wild_filepath,width = 8,height = 6)
gseaplot2(res,
          geneSetID = 1,
          title = "Splice vs Wild HALLMARK_TP53_PATHWAY",
          pvalue_table = T
)
dev.off()

#splice vs coding
res <- results(dds,contrast = c('tp53_type', "splice","coding"))
write.table(res,deseq2_tp53splice_vs_coding.txt_filepath,sep="\t", quote=F, col.names=NA) #保留中间结果
#splice vs coding-GSEA plot
tp53_geneSets <- read.gmt(hallmark_tp53_pathway_filepath)
df <- read.delim(deseq2_tp53splice_vs_coding.txt_filepath,row.names = 1)
df <- na.omit(df)
df_sort <- df %>%
  mutate(symbol = sapply(rownames(df), function(x) strsplit(x, "\\|")[[1]][2])) %>%
  arrange(desc(log2FoldChange))
gene_list <- df_sort$log2FoldChange
names(gene_list) <- df_sort$symbol
gene_list <- gene_list[gene_list!=0]
gene_list <- gene_list[!duplicated(names(gene_list))]
res <- GSEA(
  gene_list,
  TERM2GENE = tp53_geneSets,
  pvalueCutoff = 1
)
pdf(file = GSEA_tp53splice_vs_coding_filepath,width = 8,height = 6)
gseaplot2(res,
          geneSetID = 1,
          title = "Splice vs Coding HALLMARK_TP53_PATHWAY",
          pvalue_table = T
)
dev.off()

##Survival analysis
singlehit_samples <- tp53_df %>%
  mutate(type = ifelse(is.na(type),"Wild-type",
                       ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"FS",
                              ifelse(type=="Missense_Mutation","MIS",
                                     ifelse(type=="Nonsense_Mutation","NS",
                                            ifelse(type=="Splice_Mutation","SP",
                                                   ifelse(type=="In_Frame_Del","IFD",type))))))) %>%
  select(-coding_or_splice) %>%
  group_by(hugo_symbol,sample) %>%
  summarise(type = paste(type, collapse = ";"),
            .groups = "drop") %>%
  filter(type %in% c("MIS","FS","NS","SP","IFD")) %>%
  select(sample)
singlehit_samples <- singlehit_samples$sample

df <- config %>%
  mutate(sample = paste0(V1,"LC")) %>%
  left_join(tp53_df,by=c("sample")) %>%
  mutate(type = ifelse(is.na(type),"Wild-type",
                       ifelse(type %in% c("Frame_Shift_Del","Frame_Shift_Ins"),"Frame_Shift",
                              ifelse(type=="Missense_Mutation","Missense",
                                     ifelse(type=="Nonsense_Mutation","Nonsense",
                                            ifelse(type=="Splice_Mutation","Splice",type)))))) %>%
  filter(sample %in% singlehit_samples | type=="Wild-type") %>% #去除有多个TP53突变的患者
  filter(type != "In_Frame_Del") %>%#去除In_Frame_Del的患者，例数太少
  left_join(clin,by=c("sample")) %>%
  select(sample,type,OS,OS_E,RFS,RFS_E)
df$type <- factor(df$type)
#os
pdf(tp53_mutationtype_os_km_filepath,width = 8,height = 5,onefile = F)
km_os<-survfit(Surv(OS,OS_E)~type,data = df)
p <- ggsurvplot(km_os, 
                conf.int=F, #是否显示生存率的95%CI
                #risk.table=T,#显示风险表
                palette=c("#99CCFF","#E5CCFF",
                          "#CCFF99","#F2B2A2",
                          "darkred"), #配色
                title=paste0("TP53_mutation_type","_OS"), #大标题
                #risk.table.height = 0.25,#风险表的高度比例
                pval=T,  
                pval.method=T
)
print(p)
dev.off()
#rfs
pdf(tp53_mutationtype_rfs_km_filepath,width = 8,height = 5,onefile = F)
km_rfs<-survfit(Surv(RFS,RFS_E)~type,data = df)
p1 <- ggsurvplot(km_rfs, 
                 conf.int=F, #是否显示生存率的95%CI
                 #risk.table=TRUE,#显示风险表
                 palette=c("#99CCFF","#E5CCFF",
                           "#CCFF99","#F2B2A2",
                           "darkred"), #配色
                 title=paste0("TP53_mutation_type","_RFS"), #大标题
                 #risk.table.height = 0.25, #风险表的高度比例
                 pval=T,  
                 pval.method=T)
print(p1)
dev.off()







#draft code
if(0){
ssgsea<- gsva(dat, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
fgseaRes <- fgsea(pathways = tp53_geneSets, 
                  stats = gene_list,
                  eps = 0.0,
                  minSize = 15,
                  maxSize = 500)
plotEnrichment(tp53_geneSets[["gene"]],
               gene_list) + labs(title="HALLMARK_TP53_PATHWAY")
#纵坐标取对数
scale_y_log10()+
  , title = res$Description[1], geneSetID = 1
  df_1 <- df %>%
  filter(cn==0)

  str_detect("TP53",rownames(expr))
  col_sums <- colSums(mat)
  sorted_col_names <- names(sort(col_sums,decreasing = T))
  mat_sorted <- mat[,sorted_col_names]
samples1 <- unique(coding_tp53_df$sample)
samples2 <- unique(splice_tp53_df$sample)
df <- data.frame(row1 = c(samples1, rep(NA, length(samples2) - length(samples1))),
                row2 = c(samples2, rep(NA, length(samples1) - length(samples2))))

rownames <- c("TP53_coding","TP53_splice")
all_sample
mat_forplot <- 
}

