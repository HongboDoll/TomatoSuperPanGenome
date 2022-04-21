#!/bin/bash

ref=S_galapagense.fa

threads=40
for i in TS-1 TS-2
do
        ./vg giraffe -t $threads -f ${i}_1.clean.fastq.gz -f ${i}_2.clean.fastq.gz -x tomato_PAV.xg -g tomato_PAV.gg -m tomato_PAV.min -d tomato_PAV.dist -H tomato_PAV.gbwt > ${i}.gam
#
        ./vg pack -t $threads -Q 5 -x  tomato_PAV.xg -g ${i}.gam -o ${i}.pack
        ./vg call -t $threads -r tomato_PAV.snarls -k ${i}.pack  -s ${i} tomato_PAV.xg  > ${i}.vcf  #### no -v $vcf parameter, drastically reduce the memory cost
done

ls vcf | while read i
do
    cat <(grep '#' vcf/$i) <(./filter_vcf_ad.py vcf/$i | grep -v '#') > vcf/${i}.filter
    bcftools norm -m -both vcf/${i}.filter | bcftools norm -d none -c x --fasta-ref $ref | bcftools sort | bgzip > vcf/${i}.filter.norm.vcf.gz
    tabix -f vcf/${i}.filter.norm.vcf.gz
done
	
VCFS=`ls vcf/*.filter.norm.vcf.gz`
bcftools merge -0 $VCFS | bcftools norm -m -any -N -c x --fasta-ref $ref > 321_acc.merge.vcf

grep -v '#' 321_acc.merge.vcf | grep -v 'HiC_scaffold' | sort -k1,1 -k2,2n > t
cat <(grep '#' 321_SV.merge.vcf | grep -v HiC_scaffold) <(paste <(cut -f '1-2' t) <(awk 'BEGIN{n=1}{print "SV"n;n+=1}' t) <(cut -f '4-999' t)) > tt  && mv tt 321_SV.merge.vcf
