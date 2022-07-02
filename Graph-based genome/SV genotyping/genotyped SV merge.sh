#!/bin/bash

ref=S_galapagense.fa

#ls vcf | while read i
#do
#    cat <(grep '#' vcf/$i) <(./filter_vcf_ad.py vcf/$i | grep -v '#') > vcf/${i}.filter
#    bcftools norm -m -both vcf/${i}.filter | bcftools norm -d none -c x --fasta-ref $ref | bcftools sort | bgzip > vcf/${i}.filter.norm.vcf.gz
#    tabix -f vcf/${i}.filter.norm.vcf.gz
#done
#	
VCFS=`ls vcf/*.filter.norm.vcf.gz`
bcftools merge $VCFS > 321_acc.merge.vcf

grep -v '#' 321_acc.merge.vcf | grep -v 'HiC_scaffold' | sort -k1,1 -k2,2n > t
cat <(grep '#' 321_SV.merge.vcf | grep -v HiC_scaffold) <(paste <(cut -f '1-2' t) <(awk 'BEGIN{n=1}{print "SV"n;n+=1}' t) <(cut -f '4-999' t)) > tt  && mv tt 321_SV.merge.vcf

cat 321_SV.merge.vcf| grep -v '#'|awk '{print $1,$2}'|sort | uniq -c |awk '$1!=1' > duplicated_pos.xls

./remove_duplicated_sites_in_vcf.py duplicated_pos.xls 321_SV.merge.vcf > 321_SV.merge_revise.vcf

#cat 12sol_ZY56_SV_syri_merge*wild.vcf  |grep -v '#'|awk 'BEGIN{OFS="\t"}{print $1,$2-25,$2+25}' | sort -k1,1 -k2,2n >  12sol_ZY56_wild_SV_INS_DEL_INV_coordinates.xls
#cat 321_SV.merge.vcf | grep --color=auto -v '#' | awk 'BEGIN{OFS="\t"}{print $1,$2-25,$2+25}' | sort -k1,1 -k2,2n > 321_SV_coordinates.xls

#cat 12sol_ZY56_SV_syri_merge*.vcf  |grep -v '#'|awk 'BEGIN{OFS="\t"}{print $1,$2-25,$2+25}' | sort -k1,1 -k2,2n >  12sol_ZY56_SV_INS_DEL_INV_coordinates.xls
#
#bedtools intersect -v -a 321_SV_coordinates.xls -b 12sol_ZY56_SV_INS_DEL_INV_coordinates.xls > additional_sv_adding_cell.xls
#
#bedtools intersect -v -a 321_SV_coordinates.xls -b 12sol_ZY56_wild_SV_INS_DEL_INV_coordinates.xls > additional_sv_adding_cell.xls2
#
#bedtools intersect -a 321_SV_coordinates.xls -b 12sol_ZY56_SV_INS_DEL_INV_coordinates.xls | wc -l
#bedtools intersect -a 321_SV_coordinates.xls -b 12sol_ZY56_wild_SV_INS_DEL_INV_coordinates.xls |wc -l
#
