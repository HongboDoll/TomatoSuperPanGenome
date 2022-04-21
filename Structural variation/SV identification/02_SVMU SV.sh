#!/bin/bash

for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65 ZY59 Spenn ZY63 ZY64
do
	grep CNV sv.${i}_ZY56_svmu.txt | awk '$(NF-2)>50' | awk '$NF!=0&&$(NF-1)!=0&&$NF!="inf"&&$(NF-1)!="inf"' | awk '$NF/$(NF-1)>=2||$(NF-1)/$NF>=2' > sv.${i}_ZY56_svmu_CNV.xls
        rm sv.${i}_ZY56_svmu_CNV.vcf
        cat sv.${i}_ZY56_svmu_CNV.xls | while read n
        do
                chr=`echo $n | awk '{print $1}'`
                pos=`echo $n | awk '{print $2}'`
		ref_pos=`echo $n | awk '{if($2<=$3){print $1":"$2"-"$3}else{print $1":"$3"-"$2}}'`
		alt_pos=`echo $n | awk '{if($6<=$7){print $5":"$6"-"$7}else{print $5":"$7"-"$6}}'`
                ref_seq=`samtools faidx ZY56.fa $ref_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                alt_seq=`samtools faidx ${i}.fa $alt_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                echo -e "${chr}\t${pos}\t.\t${ref_seq}\t${alt_seq}\t30\tPASS\t.\tGT\t1/1" >> sv.${i}_ZY56_svmu_CNV.vcf
        done
	./uniq_vcf_chr_pos.py sv.${i}_ZY56_svmu_CNV.vcf > o; mv o sv.${i}_ZY56_svmu_CNV.vcf

done

