#!/bin/bash

for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65 ZY59 Spenn ZY63 ZY64
do
	for sv in CNV 
	do
	grep -v 'N' sv.${i}_ZY56_svmu_${sv}.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,"DP=100","'"$i"'"}' > ${i}_ZY56_${sv}.vcf1
	done
done

##### caution: vcf for each accession must be sorted and uniq ####
for sv in CNV 
do
        rm 12sol_ZY56_${sv}.vcf
        cat *_ZY56_${sv}.vcf1 > 12sol_ZY56_${sv}.vcf
        ./merge_vcf.py 12sol_ZY56_${sv}.vcf > 12sol_ZY56_${sv}.merge.vcf
        sort -k1,1 -k2,2n 12sol_ZY56_${sv}.merge.vcf -o 12sol_ZY56_${sv}.merge.vcf
done

for i in CNV 
do
cat 12sol_ZY56_${i}.merge.vcf|awk '{print $1"~"$2}'|sort |uniq -c|awk '$1!=1'|awk '{print $2}'|awk -F '~' '{print $1,$2}' > 12sol_ZY56_${i}.problem.pos.xls
done

### some indels have the identical reference coordinates but different allels, which should be addressed

for i in CNV 
do
split -a 2 -d -l 140 12sol_ZY56_${i}.problem.pos.xls 12sol_ZY56_${i}.problem.pos.xls_
done


for n in {00..24}
do
for i in CNV 
do
rm 12sol_ZY56_${i}.problem.pos.correct.vcf
echo -e "cat 12sol_ZY56_${i}.problem.pos.xls_${n} | while read i;     do         echo \$i > tmp_${i}_${n};         ./correct_identical_reference_pos_indel.py 12sol_ZY56_${i}.merge.vcf tmp_${i}_${n} >> 12sol_ZY56_${i}.problem.pos.correct.vcf;     done" > correct_indel.sh_${i}_${n}
chmod 755 correct_indel.sh_${i}_${n}
done
done

for n in {00..24}
do

for i in CNV 
do
qsub -N corr_${i}_${n} -pe smp 1 -S /bin/bash -cwd -j y -V -b y "./correct_indel.sh_${i}_${n}"
done
done

###############################
###############################
for i in CNV 
do
./remove_problematic_indel_in_vcf.py 12sol_ZY56_${i}.problem.pos.xls 12sol_ZY56_${i}.merge.vcf > tmp
cat 12sol_ZY56_${i}.problem.pos.correct.vcf tmp | sort -k1,1 -k2,2n > tt && mv tt 12sol_ZY56_${i}.merge.vcf
done

cat *_svmu_CNV.xls > tmp
./add_cnv_id.py tmp 12sol_ZY56_CNV.merge.vcf > tmp1; mv tmp1 12sol_ZY56_CNV.merge.vcf

