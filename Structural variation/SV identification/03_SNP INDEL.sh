#!/bin/bash

nucmer -t 40 $ref3 $ref2 --prefix ZY64_ZY56
delta-filter -1 ZY64_ZY56.delta > ZY64_ZY56.filter.delta

ls *.delta | while read i
do
i=ZY64_ZY56.delta
query=`echo $i | awk -F '.' '{print $1}' | awk -F '_' '{print $1}'`
lastz ZY56.fa ${query}.fa --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > ${query}_ZY56.sam_lastz.txt
qsub -N ${query} -l p=1 -l h_vmem=80g  -q all.q -S /bin/bash -cwd -j y -V -b y  "svmu $i ZY56.fa ${query}.fa h ${query}_ZY56.sam_lastz.txt ${query}_ZY56_svmu"
done

ls *_ZY56.delta | while read i
do
i=ZY64_ZY56.delta
nucmer -t 20 ZY56.fa ZY64.fa -p ZY64_ZY56
uery=`echo $i | awk -F '.' '{print $1}' | awk -F '_' '{print $1}'`
svmu $i ZY56.fa ${query}.fa h ${query}_ZY56.sam_lastz.txt ${query}_ZY56_svmu
./filter_delta_based_svmu_cm.py cm.${query}_ZY56_svmu.txt $i > ${query}_ZY56_svmu_coords_filter.delta
done

ls *_ZY56.delta | while read i
do
query=`echo $i | awk -F '.' '{print $1}' | awk -F '_' '{print $1}'`
show-snps -Clr -x 1 -T ${query}_ZY56_svmu_coords_filter.delta > ${query}_ZY56.snv
./MUMmerSNPs2VCF.py ${query}_ZY56.snv ${query}_ZY56_snv.vcf
all=`show-coords -THrcl ${query}_ZY56_svmu_coords_filter.delta |awk '{i+=$7}END{printf "%.f",i}'`
num=`show-coords -THrcl ${query}_ZY56_svmu_coords_filter.delta | wc -l`
ave_iden=`expr $all \/ $num `
echo $query $ave_iden
done

