#!/bin/bash 
#===============================================================================
#          FILE: VEPoutput(jiaoji) to enpm_adw_mutspot input
#         USAGE: ./VEPO_jiaoji_2_enpm_adw_mutspot_input.sh
#   DESCRIPTION: 这个脚本用于去除VEP注释文件中的SNP，然后得到enpm,adw,mutspot三个软件所需的input
#       OPTIONS: ---
#  REQUIREMENTS: 
#          BUGS: ---
#         NOTES: 主要步骤是去snp，去重复行(vep一个突变有多条注释);另外适用于hg38，如果是hg19需要liftover
#        AUTHOR: 唐华韬
#  ORGANIZATION: 
#       CREATED: 1/14/25
#      REVISION:  ---
#      Reference:
#===============================================================================
set -o nounset                              # Treat unset variables as an error

vepoutput_path=~/vep_data/986_jiaoji.vepoutput
noncoding_snv_analysis_path=/home/tht/noncoding_1000pair/3.ssm_driver_calling/noncoding_snv_analysis
remove_snp=/home/tht/noncoding_1000pair/3.ssm_driver_calling/noncoding_snv_analysis/remove_snp.R
activedriverwgs_path=/home/tht/noncoding_1000pair/3.ssm_driver_calling/activedriverwgs
mutspot_path=/home/tht/noncoding_1000pair/3.ssm_driver_calling/mutspot/MutSpot/
enpminput_filename="986_jiaoji.enpm.input"
countvsrecn_filename="986_jiaoji.cr.input"
adwinput_filename="986_jiaoji.adwmuts.input"
mutspot_snv_filename="986_jiaoji.mutspot_snv.input"
mutspot_indel_filename="986_jiaoji.mutspot_indel.input"

#去掉所有##的行
cd $noncoding_snv_analysis_path
cp $vepoutput_path vepoutput
grep -v '^##' vepoutput > output_file

#1 得到enpm和countvsrecn的input
#去掉snp 
Rscript $remove_snp $noncoding_snv_analysis_path/output_file
#取出需要的列并去除重复的行
cat vep_result_germline_removed.txt | awk 'NR!=1{split($1, a, ":");OFS="\t";print a[1],a[2],a[3],a[4],a[5]}' | sort | uniq > $countvsrecn_filename
rm -f vepoutput output_file

#取出重复次数大于等于n次的行
awk  '{key = $2 OFS $3; count[key]++; lines[key] = lines[key] $0 ORS} END {for (k in count) if (count[k] >= 10 && lines[k] != "") print lines[k]}' OFS='\t' $countvsrecn_filename> enpm.input1
#删除空行
awk 'NF' enpm.input1 > $enpminput_filename

#去除中间文件
rm -f vep_result_germline_removed.txt enpm.input1 

#2 得到ADWinput
#这个软件的输入表格格式是chr,突变起始位置(1-based)，突变结束位置(1-based)，ref，alt，patient
cd $activedriverwgs_path
cp $noncoding_snv_analysis_path/$countvsrecn_filename crinput
./crinput2ADWinput.awk crinput > $adwinput_filename

rm -f crinput

#得到 MutSpot input
#这个软件的输入表格格式和ADW一样，只是snv必须和indel分开
awk 'NR!=1 && $4 ~ /^[a-zA-Z]$/ && $5 ~ /^[a-zA-Z]$/ {OFS="\t";print $0}' $adwinput_filename > $mutspot_path/$mutspot_snv_filename
awk 'NR!=1 && ($4 !~ /^[a-zA-Z]$/ || $5 !~ /^[a-zA-Z]$/) {OFS="\t";print $0}' $adwinput_filename> $mutspot_path/$mutspot_indel_filename



