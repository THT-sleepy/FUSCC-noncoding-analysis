#!/bin/bash 
#===============================================================================
#          FILE: vcftovepoutput_for986fdscc_wgs_pairs_all.sh
# 
#         USAGE: sudo ./vcftovepoutput_for986fdscc_wgs_pairs_all.sh
# 
#   DESCRIPTION: 这个脚本用于把986个样本的VCF用于VEP注释，得到VEPoutput,对VCF的筛选是所有,但这个脚本31行开始需要手动run
# 
#       OPTIONS: ---
#  REQUIREMENTS: VEP的docker;run.sh(在脚本末尾);vcf文件;vcf2vepinput_all.sh(在脚本末尾);config文件
#          BUGS: ---
#         NOTES: 
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
export VCF2VEPINPUT=/mnt/sdc/tanght/noncoding_1000pair/2.annotation_vep/vcf2vepinput_all.sh
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
wait
#合并起来
grep -v "^#" $vepinput_part2_filename.vepoutput > ${vepinput_part2}_noheader
grep -v "^#" $vepinput_part3_filename.vepoutput > ${vepinput_part3}_noheader
cat $vepinput_part1_filename.vepoutput ${vepinput_part2}_noheader ${vepinput_part3}_noheader> $vepoutput

#删除中间文件
rm -f $vepinput_part1_filename.vepoutput ${vepinput_part2}_noheader ${vepinput_part3}_noheader  $vepinput_part2_filename.vepoutput $vepinput_part3_filename.vepoutput  $vepinput_part1_filename $vepinput_part2_filename $vepinput_part3_filename 

#把986_all.vepinput移动回来
mv $vepinput_filename $annotation_vep_dir

#run.sh
<<'COMMENT'
time sudo docker run -v /mnt/sdc/tanght/vep_data:/data ensemblorg/ensembl-vep \
        vep --cache \
        --input_file /data/$1 \
        --output_file /data/$1.vep \
        --tab \
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
        --gencode_basic
COMMENT

#vcf2vepinput.sh
<<'COMMENT'
#vcf2vepinput.sh
#get chr1-22,X,Y mutations supported by more than one ssm caller
zcat $VCF/$1*.vcf.gz | perl -alne '{next if $F[0] =~ /_/ || $F[0] =~ /^#/ ; print unless $F[0] =~ /V/ }' > $1.vcf1
cat $1.vcf1 | awk -v id=$1 '{OFS="\t";var=id ":" $1 ":" $2 ":" $4 ":" $5 ":" $8;print $1,$2,var,$4,$5,$6,$7,$8}' >> 986.>
rm $1.vcf1
COMMENT



