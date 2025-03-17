#vcf2vepinput.sh
#get chr1-22,X,Y mutations supported by more than one ssm caller
zcat $VCF/$1*.vcf.gz | perl -alne '{next if $F[0] =~ /_/ || $F[0] =~ /^#/ ; print unless $F[0] =~ /V/ }' > $1.vcf1
cat $1.vcf1 | awk -v id=$1 '{OFS="\t";var=id ":" $1 ":" $2 ":" $4 ":" $5 ":" $8;print $1,$2,var,$4,$5,$6,$7,$8}' >> 986.merge.vcf1 #set id and add format column
rm $1.vcf1
