# Quality check on FUSCC SV data             (What)

* Aug 23, 2025                                 (When)
* biotrainee ~/1000_noncoding/3.SV_analysis/2.analysis           (Where)
* 对1000对的SV输出做一个质控(PR>=3,SR>=1)         (Why)


## 代码                                          (How)

### loading packages
```
library(vcfR)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(purrr)
```

### set variables
```
folder_path <- "~/1000_noncoding/3.SV_analysis/1.resources/manta_results/somaticSV_converInved/"
file_list <- list.files(path = folder_path, pattern = "\\.gz$", full.names = TRUE)
purity_path <- "~/1000_noncoding/信息源文件/tumor_purity_20230405.csv"
```
### get merged vcf
```
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
all_vcf_data <- lapply(file_list, process_vcf) %>%
  discard( ~ nrow(.x) == 0) %>%
  bind_rows()
save(all_vcf_data,file="~/1000_noncoding/3.SV_analysis/2.analysis/RData/lc986_vcf_afterQC.Rdata")
```


### 其它文章LUAD sv每个样本看了下大概100出头这样子，我们的PR,SR标准严格一点均值是54（Result）

# to do list
## 0 要修正一下sv的计数(reciprocal算两个？)-done
## 1 比较一下不同分期总sv数量的差别以及各类型sv的差别(箱线图+百分比柱状图)-done
## 2 sv 和 cnv的一致性(cnv要用GISTC2对应的那个结果)-done
## 3 sv类型和复制时间早晚的关系-done
## 4 sv数量和临床指标的关系，和驱动突变以及TP53等突变的关系-done
## 5 sv数量类型和染色质可及性的关系-done
## 6 搞一个多重回归，看下sv数量到底受什么影响比较大
