# Quality check on FUSCC SV data             (What)

* Aug 23, 2025                                    (When)
* biotrainee ~/1000_noncoding/3.SV_analysis/2.analysis           (Where)
* 对1000对的SV输出做一个质控,主要是PR和SR                (Why)


## 代码                                          (How)

```
##Quality control of SV variants of 1000 cohort##

##loading packages
library(vcfR)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(purrr)

##set variables
folder_path <- "~/1000_noncoding/3.SV_analysis/1.resources/manta_results/somaticSV_converInved/"
file_list <- list.files(path = folder_path, pattern = "\\.gz$", full.names = TRUE)
purity_path <- "~/1000_noncoding/信息源文件/tumor_purity_20230405.csv"

##1 get merged vcf
#定义操作每个vcf文件的函数
process_vcf <- function(file) {
  # 使用 read.vcfR 读取 VCF 文件，得到 S4 对象
  vcf <- read.vcfR(file)

  # 提取 fix 部分并转换为数据框
  gt <- vcf@gt
  vcf <- vcf@fix
  vcf <- as.data.frame(vcf)
  gt <- as.data.frame(gt)
  vcf <- cbind(vcf,gt)
  sample <-sub("(.*LC).*", "\\1", colnames(gt)[3])
  vcf <- vcf %>%
    #去除PR<3的以及FILTER == "PASS" 的记录
    mutate(PR=sapply(vcf[,c(11)], function(x) as.numeric(strsplit(x, "[,:]")[[1]][2])))%>%
    filter(PR>3,FILTER == "PASS") %>%
    #去除SR为0的
    filter(FORMAT != "PR") %>%
    mutate(
      SV_ID = sapply(str_split(ID, ":"), function(x) {
        if ("TANDEM" %in% x) {
          # 如果包含TANDEM，取前三个值，并在前面加上 sample
          paste(sample, paste(x[1:3], collapse = ":"), sep = ":")
        } else {
          # 否则取前两个值，并在前面加上sample
          paste(sample, paste(x[1:2], collapse = ":"), sep = ":")
        }
      }),
      SV_TYPE = sapply(strsplit(INFO, ";"), function(x) {
        # 提取SVTYPE的值
        sv_type_value <- sub("SVTYPE=([A-Za-z]+)", "\\1", x[grepl("SVTYPE", x)])
        # 根据条件修改SV_TYPE
        if (sv_type_value == "BND") {
          return("TRA")
        } else {
          return(sv_type_value)
        }
      }),SAMPLE=sample
    )
vcf <- vcf[,-c(9,10,11,12)]
  return(vcf)
}

#使用 lapply 读取所有文件并合并成一个数据框
load("~/1000_noncoding/3.SV_analysis/2.analysis/lc986_vcf_afterQC.Rdata")
if(0){
all_vcf_data <- lapply(file_list, process_vcf) %>%
  discard( ~ nrow(.x) == 0) %>%
  bind_rows()
}
##2 画一个排序柱状图
df <- all_vcf_data[,c("SAMPLE","SV_TYPE")]
df$SV_TYPE <- factor(df$SV_TYPE, levels = c("TRA","INS","INV","DUP","DEL"))
p1_1 <- ggplot(df, aes(x= reorder(SAMPLE, SAMPLE, function(x) -length(x)),fill = SV_TYPE)) +
  geom_bar() +
  labs(x="Tumors",y = "SomaticSVs",tag = "A") +
  theme_minimal()+
  scale_fill_manual(values = c(
    "DEL" = rgb(226, 145, 53, maxColorValue = 255),
    "DUP" = rgb(114, 176, 99, maxColorValue = 255),
    "INV" = rgb(113, 154, 172, maxColorValue = 255),
    "INS" = rgb(184, 219, 179, maxColorValue = 255),
    "TRA" = rgb(74, 95, 126, maxColorValue = 255)
  )) +
  theme(axis.text.x = element_blank(),
        plot.tag = element_text(size = 16, face = "bold"), #标签加粗并增大字体
        axis.title.x = element_text(size = 16, face = "bold"), # 横坐标标题字体加大加粗
        axis.title.y = element_text(size = 16, face = "bold"), # 纵坐标标题字体加大加粗
        legend.position = c(0.8, 0.6))

##3 查看各个样本的count,画一个各类型的柱状图
df <- all_vcf_data %>%
  group_by(SAMPLE) %>%
  count() %>%
  arrange(desc(n))
summary(df)
#SAMPLE                n
#Length:894         Min.   :   1.00
#Class :character   1st Qu.:  13.00
#Mode  :character   Median :  40.00
#Mean   :  76.09
#3rd Qu.: 102.00
#Max.   :2645.00
df <- all_vcf_data[,c("SAMPLE","SV_TYPE")]
p1_2 <- ggplot(df, aes(x = SV_TYPE,fill = SV_TYPE)) +
  geom_bar() +
  scale_fill_manual(values = c(
    "DEL" = rgb(226, 145, 53, maxColorValue = 255),
    "DUP" = rgb(114, 176, 99, maxColorValue = 255),
    "INV" = rgb(113, 154, 172, maxColorValue = 255),
    "INS" = rgb(184, 219, 179, maxColorValue = 255),
    "TRA" = rgb(74, 95, 126, maxColorValue = 255)
  ))+
  labs(tag = "B",x="",y="",title = "") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, color = "black") + # 添加文本标签
  theme_minimal()+
  ylim(c(0,31000))+
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),  # 标题加粗并增大字体
    axis.title.x = element_text(size = 14, face = "bold"),  # 横坐标加粗并增大字体
    axis.title.y = element_text(size = 14, face = "bold"),   # 纵坐标加粗并增大字体
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标字体加大加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标字体加大加粗
    legend.position = "none"
  )


##4 查看purity和count的关系-相关性曲线
purity <- read.csv(purity_path)
df <- all_vcf_data %>%
  group_by(SAMPLE) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(sample_id=SAMPLE) %>%
  left_join(purity,by =c("sample_id"))
df1 <- na.omit(df)

# 计算相关系数和p值
cor_test <- cor.test(df1$purity, df1$n)
correlation <- cor_test$estimate  # 相关系数
p_value <- cor_test$p.value      # p值

p1_3 <- ggplot(df1, aes(x = purity, y = n)) +
  geom_point(color = "darkblue") +  # 绘制散点
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # 拟合直线
  annotate("text", x = 0.8, y = 3000,
           label = paste("Correlation: ", round(correlation, 2),
                         "\np-value: ", format(p_value, digits = 3)),
           color = "black", size = 5) +  # 显示相关系数p
  labs(tag="C",title = "Correlation between purity and somaticSVs", x = "purity", y = "somaticSVs") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题加粗并增大字体
    plot.tag = element_text(size = 16, face = "bold"), #标签加粗并增大字体
    axis.title.x = element_text(size = 14, face = "bold"),  # 横坐标加粗并增大字体
    axis.title.y = element_text(size = 14, face = "bold"),   # 纵坐标加粗并增大字体
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标字体加大加粗
    axis.text.y = element_text(size = 12, face = "bold")  # 纵坐标字体加大加粗
  )


##5 去掉hypermutated以及肿瘤纯度过低的样本并保存对象
#先看一下突变数量的分布情况
purity <- read.csv(purity_path)
df <- all_vcf_data %>%
  group_by(SAMPLE) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(sample_id=SAMPLE) %>%
  left_join(purity,by=c("sample_id"))

p1_4 <- ggplot(df, aes(x = n)) +
  geom_histogram(binwidth = 30,  # 控制区间宽度
                 fill = "darkred" # 填充颜色
                 ) +
  labs(title = "",      # 标题
       x = "SomaticSVs",                 # x 轴标签
       y = "Count",
       tag = "D")+
  theme_minimal()+
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),  # 标题加粗并增大字体
    axis.title.x = element_text(size = 14, face = "bold"),  # 横坐标加粗并增大字体
    axis.title.y = element_text(size = 14, face = "bold"),   # 纵坐标加粗并增大字体
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标字体加大加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标字体加大加粗
  )
# 计算四分位数和 IQR
Q1 <- quantile(df$n, 0.25)  # 下四分位数
Q3 <- quantile(df$n, 0.75)  # 上四分位数
IQR <- Q3 - Q1  # 四分位距

# 计算上下限
lower_limit <- Q1 - 1.5 * IQR
upper_limit <- Q3 + 1.5 * IQR

df1 <- na.omit(df)
df1 <- df1[df1$n < upper_limit,]
p1_5 <- ggplot(df1, aes(x = n)) +
  geom_histogram(binwidth = 5,  # 控制区间宽度
                 fill = "darkred" # 填充颜色
  ) +
  labs(title = "",      # 标题
       x = "SomaticSVs",                 # x 轴标签
       y = "Count",
       tag = "E")+
  theme_minimal()+
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),  # 标题加粗并增大字体
    axis.title.x = element_text(size = 14, face = "bold"),  # 横坐标加粗并增大字体
    axis.title.y = element_text(size = 14, face = "bold"),   # 纵坐标加粗并增大字体
    axis.text.x = element_text(size = 12, face = "bold"),  # 横坐标字体加大加粗
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵坐标字体加大加粗
  )
all_vcf_data_filtered <- all_vcf_data[all_vcf_data$SAMPLE %in% df1$SAMPLE,]
save(all_vcf_data,file = "~/1000_noncoding/3.SV_analysis/2.analysis/lc986_vcf_afterQC.Rdata")
save(all_vcf_data_filtered,file = "~/1000_noncoding/3.SV_analysis/2.analysis/all_vcf_afterQC.Rdata")

##6 简简单单用patchwork拼个图，复杂一定的成品图还是拼接还得是ai
p1_1 / (p1_2 | p1_3) / (p1_4 | p1_5)

```
                                                    (Result)
