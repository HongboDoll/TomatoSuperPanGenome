#!/bin/bash

plink --vcf 321_SV.merge.vcf --recode12 --allow-extra-chr --geno 0.1 --mac 5 --allow-no-sex --biallelic-only --out 321_tomato_SV

cat 321_tomato_SV.map | sed 's/Sly0//g' | sed 's/Sly//g' > 321_tomato_SV.map2
mv 321_tomato_SV.map2 321_tomato_SV.map

# PCA
plink --allow-extra-chr --file 321_tomato_SV --indep-pairwise 50 5 0.2 --recode vcf-iid --out 321_tomato_SV_LDpruned

plink --allow-extra-chr --file 321_tomato_SV --recode vcf-iid --extract 321_tomato_SV_LDpruned.prune.in --out 321_tomato_SV_LDpruned

plink --allow-extra-chr --vcf 321_tomato_SV_LDpruned.vcf --make-bed --out 321_tomato_SV_LDpruned_bfile

plink --bfile 321_tomato_SV_LDpruned_bfile --allow-extra-chr --allow-no-sex --pca 10 --out PCA

cat PCA.eigenvec | awk '{print $1,$2,"1",$3,$4,$5,$6,$7}' > PCA_SV.txt

threshold=`grep -v "Observed_Number" 321_tomato_SV.gec.sum | awk '{print $4}'`

# EMMAX
plink --file 321_tomato_SV --allow-extra-chr --recode 12 --output-missing-genotype 0 --transpose --out 321_tomato_SV

emmax-kin -v -d 10 321_tomato_SV

rm -rf sv_results; mkdir sv_results
ls phenotype/*xls | while read i
do
trait=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
emmax -v -d 10 -t 321_tomato_SV \
        -p $i \
        -k 321_tomato_SV.aBN.kinf \
        -c PCA_SV.txt \
        -o sv_results/$trait

cat 321_tomato_SV.map | paste - sv_results/${trait}.ps | awk '{print $2,$1,$4,$8}' | sed '1i\SNP CHR BP P' | sed 's/ /\t/g' > sv_results/${trait}.321_tomato_SV.out.txt

./manhattan_qq.R sv_results/${trait}.321_tomato_SV.out.txt sv_results/${trait}_EMMAX_aBN.ps.png $threshold
done

