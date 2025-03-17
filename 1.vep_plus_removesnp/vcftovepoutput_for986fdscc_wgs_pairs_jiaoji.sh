#!/bin/bash 
#===============================================================================
#          FILE: vcftovepoutput_for986fdscc_wgs_pairs_jiaoji.sh
# 
#         USAGE: sudo ./vcftovepoutput_for986fdscc_wgs_pairs_jiaoji.sh
# 
#   DESCRIPTION: 这个脚本用于把986个样本的VCF用于VEP注释，得到VEPoutput
#       OPTIONS: ---
#  REQUIREMENTS: VEP的docker;run.sh(在脚本末尾)
#          BUGS: ---
#         NOTES: 环境变量里面文件路径整复杂了，下次应该把文件夹和文件名分开设置会清楚一些，另外可以添加一些执行到哪了的echo 命令
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 7/12/24 
#      REVISION:  ---
#      Reference: 
#===============================================================================
set -o nounset                              # Treat unset variables as an error
#variables
export config=/mnt/sdc/tanght/noncoding_1000pair/1.resources/986.config
export VCF=/mnt/sdc/tanght/noncoding_1000pair/1.resources/986vcfs
export VCF2VEPINPUT=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/vcf2vepinput_jiaoji.sh
export RUN_VEP=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/run.sh
export merge_vcf1=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986.merge.vcf1
export merge_vcf2=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986.merge.vcf2
export vepinput=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput
export vepinput_filename=986_all.vepinput
export annotation_vep_dir=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep
export vepinput_split1=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_part_aa
export vepinput_split2=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_part_ab
export vepinput_split3=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_part_ac
export vepinput_part1=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_part1
export vepinput_part2=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_part2
export vepinput_part3=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_part3
export vepinput_part1_filename=986_all.vepinput_part1
export vepinput_part2_filename=986_all.vepinput_part2
export vepinput_part3_filename=986_all.vepinput_part3
export header=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepinput_header
export vep_filepath=/mnt/sdc/tanght/vep_data
export vepoutput=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/986_all.vepoutput
export splitpref=986_all.vepinput_part_
export TMPDIR=/mnt/sdc/tanght

#get vepinput
while read id; do $VCF2VEPINPUT "$id"  ; done < $config #get merge_vcf1 3min
cat $merge_vcf1 | awk 'BEGIN{print"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}{print $0}' > $merge_vcf2 #add header
{ grep "^#" $merge_vcf2; grep -v "^#" $merge_vcf2  | sort -k1,1 -k2,2n; } > $vepinput #sort 1-2min
rm -f $merge_vcf1 $merge_vcf2 #删除中间文件

#为了最大化利用CPU以节约时间，可以将输入三等分run
lines=$(wc -l < $vepinput) # 获取文件总行数
lines_per_file=$(((lines+2) / 3)) # 计算每个文件应包含的行数
split -l $lines_per_file $vepinput $splitpref # 使用 split 分割文件
#后面两个文件需要加上header
cat $vepinput_split1 | awk 'NR==1{print $0}' > $header 
cat $header $vepinput_split2 > $vepinput_part2 #add header
cat $header $vepinput_split3 > $vepinput_part3 #add header
mv $vepinput_split1 $vepinput_part1 #修改split1的名字为part1
rm -f $vepinput_split2 $vepinput_split3 
mv $vepinput_part1 $vepinput_part2 $vepinput_part3 $vepinput $vep_filepath #把需要的文件移动到vep文件夹下

#然后$RUN_VEP 每一部分
cd $vep_filepath
$RUN_VEP $vepinput_part1_filename > $vepinput_part1.log.txt 2>&1 &
$RUN_VEP $vepinput_part2_filename > $vepinput_part2.log.txt 2>&1 &
$RUN_VEP $vepinput_part3_filename > $vepinput_part3.log.txt 2>&1 &
wait #确保后台进程执行完再执行下一步

#合并起来
grep -v "^#" $vepinput_part2_filename.vepoutput > ${vepinput_part2}_noheader
grep -v "^#" $vepinput_part3_filename.vepoutput > ${vepinput_part3}_noheader
cat $vepinput_part1_filename.vepoutput ${vepinput_part2}_noheader ${vepinput_part3}_noheader> $vepoutput

#删除中间文件
rm -f $vepinput_part1_filename.vepoutput ${vepinput_part2}_noheader ${vepinput_part3}_noheader  $vepinput_part2_filename.vepoutput $vepinput_part3_filename.vepoutput $vepinput_part1_filename $vepinput_part2_filename $vepinput_part3_filename 

#把986_all.vepinput移动回来
mv $vepinput_filename $annotation_vep_dir

#最后把986_all.vepinput的INFO列添加到986_all.vepoutput里
grep "^#" $vepoutput > $vepoutput.title #把#开头的行先拿出来
grep "^#" $vepoutput > ${vepoutput}.1 #把不以#开头的行拿出来
grep "^#" $annotation_vep_dir/$vepinput_filename | cut -f 3,8 > $annotation_vep_dir/$vepinput_filename.1 #把vepinput不以#开头的行拿出来
sort ${vepoutput}.1 -o ${vepoutput}.1
sort $annotation_vep_dir/$vepinput_filename.1 -o $annotation_vep_dir/$vepinput_filename.1
join -a1 -1 1 -2 1 ${vepoutput}.1 $annotation_vep_dir/$vepinput_filename.1 | sed 's/\s*$/ NA/' > ${vepoutput}_withinfo1
cat $vepoutput.title ${vepoutput}_withinfo1 > ${vepoutput}_withinfo
#-a1：这个选项表示输出左表（file1.txt）中的所有行，即使在右表（file2.txt）中没有匹配。
#-1 1：表示使用左表的第一列作为连接字段（默认是第一列）。
#-2 1：表示使用右表的第一列作为连接字段（默认也是第一列）。



: <<'EOF'
#run.sh
time sudo docker run -v /mnt/sdc/tanght/vep_data:/data ensemblorg/ensembl-vep \
        vep --cache \
        --input_file /data/$1 \
        --output_file /data/$1.vepoutput \
		--format vcf \
        --variant_class \
        --sift b \
        --polyphen b \
        --regulatory \
        --cell_type A549 \
        --biotype \
        --numbers \
        --mirna  \
        --symbol \
        --tsl \
        --canonical \
        --check_existing \
        --max_af  \
        --filter_common \
        --gencode_basic \
        --buffer_size 500000 \
		--tab
#vcf2vepinput_jiaoji.sh
#get chr1-22,X,Y mutations supported by more than one ssm caller
zcat $VCF/$1*.vcf.gz | perl -alne '{next if $F[0] =~ /_/ || $F[0] =~ /^#/ ; print unless $F[0] =~ /V/ }' > $1.vcf1
cat $1.vcf1 | awk -v id=$1 '$7=="."{OFS="\t";var=id":""var"NR;print $1,$2,var,$4,$5,$6,$7,$8}' >> 986.merge.vcf1 #set id and add format column
rm $1.vcf1
#AddInfo_col.R
library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
vepinput <-fread(Args[1])
vepoutput<-fread(Args[2])
names(vepoutput)[1] <- "ID"
vepoutput<- vepoutput %>%
        left_join(vepinput,by="ID")
write.table(vepoutput, file = "986_all_canonical.vepoutput_infoadded", sep = "\t", row.names = FALSE, quote = F)
EOF
