#!/bin/bash

cat vg_vcf_head.xls <(cat <(cat 12_100_merged_ZY56_DEL_revised.vcf |grep -v '#'|awk 'BEGIN{OFS="\t";n=1}{print $1,$2,"SV_DEL"n,$4,$5,"100","PASS","DP=1000;VT=INDEL";n+=1}') <(cat 12_100_merged_ZY56_INS_revised.vcf |grep -v '#'|awk 'BEGIN{OFS="\t";n=1}{print $1,$2,"SV_INS"n,$4,$5,"100","PASS","DP=1000;VT=INDEL";n+=1}') <(cat 12_100_merged_ZY56_INV_revised.vcf |grep -v '#'|grep SV_INV|awk 'length($4)<1000000' | awk 'BEGIN{OFS="\t";n=1}{print $1,$2,"SV_INV"n,$4,$5,"100","PASS","DP=1000;VT=INDEL";n+=1}') | sort -k1,1 -k2,2n) > tomato_PAV_format.vcf
##
vcf=tomato_PAV_format.vcf
#
bgzip -f $vcf
#
tabix -f -p vcf $vcf\.gz

vcf=tomato_PAV_format.vcf.gz
ref=S_galapagense.fa
threads=40

./vg construct -t $threads -a -f -S -v $vcf -r $ref  > tomato_PAV.vg

rm -rf tmp; mkdir tmp
./vg index -t $threads -L -b tmp -x tomato_PAV.xg tomato_PAV.vg

./vg gbwt -d tmp -g tomato_PAV.gg -x tomato_PAV.xg -o tomato_PAV.gbwt -P

./vg snarls -t $threads --include-trivial tomato_PAV.xg > tomato_PAV.trivial.snarls

./vg index -b tmp -t $threads -j tomato_PAV.dist -s tomato_PAV.trivial.snarls tomato_PAV.vg

./vg minimizer -t $threads -i tomato_PAV.min -g tomato_PAV.gbwt -d tomato_PAV.dist tomato_PAV.xg

./vg snarls -t $threads tomato_PAV.xg > tomato_PAV.snarls

