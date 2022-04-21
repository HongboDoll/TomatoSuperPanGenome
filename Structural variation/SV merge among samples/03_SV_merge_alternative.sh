#!/bin/bash

ls vcf | while read i
do
    cat <(grep '#' vcf/$i) <(./filter_vcf_ad.py vcf/$i | grep -v '#') > vcf/${i}.filter
    bcftools norm -m -both vcf/${i}.filter | bcftools norm -d none -c x --fasta-ref $ref | bcftools sort | bgzip > vcf/${i}.filter.norm.vcf.gz
    tabix -f vcf/${i}.filter.norm.vcf.gz
done

VCFS=`ls vcf/*.filter.norm.vcf.gz`
sv=INS
bcftools merge -0 $VCFS | bcftools norm -m -any -N -c x --fasta-ref $ref > 12sol_ZY56_SV_syri_merge_${sv}.vcf

cat 12_sol_head <(./convert_12_sol_100_cell_vcf.py 12sol_ZY56_SV_syri_merge_${sv}.vcf) > t; mv t 12sol_ZY56_SV_syri_merge_${sv}.vcf

ls 12sol_ZY56_SV_syri_merge_${sv}.vcf > vcf_4_merge
SURVIVOR merge vcf_4_merge 50 1 0 0 0 0 merged.vcf
./revise_vcf_GT1.py merged.vcf 12sol_ZY56_SV_syri_merge_${sv}.vcf > t && mv t 12sol_ZY56_SV_syri_merge_${sv}.vcf 
./revise_merge_vcf_CHROM.py 12sol_ZY56_SV_syri_merge_${sv}.vcf head_12 > t && mv t 12sol_ZY56_SV_syri_merge_${sv}.vcf

