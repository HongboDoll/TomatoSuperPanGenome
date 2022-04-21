#!/bin/bash

for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65 ZY59 Spenn ZY63 ZY64
do
        for sv in TRANS INV  INS DEL
        do
        awk '$4!~/N/&&$5!~/N/' ${i}_ZY56_syri_${sv}.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2,".",$4,$5,"PASS",".","DP=100","'"$i"'"}' > ${i}_ZY56_${sv}.vcf1
        done
done
#
for sv in TRANS INV  INS DEL
do
        rm 12sol_ZY56_${sv}.vcf
        cat *_ZY56_${sv}.vcf1 > 12sol_ZY56_${sv}.vcf
        ./merge_vcf.py 12sol_ZY56_${sv}.vcf > 12sol_ZY56_${sv}.merge.vcf
        sort -k1,1 -k2,2n 12sol_ZY56_${sv}.merge.vcf -o 12sol_ZY56_${sv}.merge.vcf
done
#
for i in TRANS INV INS DEL
do
cat 12sol_ZY56_${i}.merge.vcf|awk '{print $1"~"$2}'|sort |uniq -c|awk '$1!=1'|awk '{print $2}'|awk -F '~' '{print $1,$2}' > 12sol_ZY56_${i}.problem.pos.xls
done
#
#### some indels have the identical reference coordinates but different allels, which should be addressed
#
for i in TRANS INV INS DEL
do
split -a 2 -d -l 440 12sol_ZY56_${i}.problem.pos.xls 12sol_ZY56_${i}.problem.pos.xls_
done

rm correct_indel*
for n in {00..02}
do
for i in TRANS DEL  INV INS
do
rm 12sol_ZY56_${i}.problem.pos.correct.vcf
echo -e "cat 12sol_ZY56_${i}.problem.pos.xls_${n} | while read i;     do         echo \$i > tmp_${i}_${n};         ./correct_identical_reference_pos_indel.py 12sol_ZY56_${i}.merge.vcf tmp_${i}_${n} >> 12sol_ZY56_${i}.problem.pos.correct.vcf;     done" > correct_indel.sh_${i}_${n}
chmod 755 correct_indel.sh_${i}_${n}
done
done

for n in {00..02}
do
for i in TRANS DEL  INV INS 
do
qsub -N corr_${i}_${n} -pe smp 1 -l h_vmem=10G -S /bin/bash -cwd -j y -V -b y "./correct_indel.sh_${i}_${n}"
./correct_indel.sh_${i}_${n} 
done
done

##############################
##############################
for i in CNV #TRANS DEL  INV INS
do
./remove_problematic_indel_in_vcf.py 12sol_ZY56_${i}.problem.pos.xls 12sol_ZY56_${i}.merge.vcf > tmp
cat 12sol_ZY56_${i}.problem.pos.correct.vcf tmp | sort -k1,1 -k2,2n > tt && mv tt 12sol_ZY56_${i}.merge.vcf
done

cat <(head -1 12sol_ZY56_DEL.merge.vcf) <(grep -v '1CHROM' 12sol_ZY56_DEL.merge.vcf |awk 'BEGIN{n=1;OFS="\t"}{print $1,$2,"SV_DEL"n,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21;n+=1}') > t1; mv t1 12sol_ZY56_DEL.merge.vcf
cat <(head -1 12sol_ZY56_INS.merge.vcf) <(grep -v '1CHROM' 12sol_ZY56_INS.merge.vcf |awk 'BEGIN{n=1;OFS="\t"}{print $1,$2,"SV_INS"n,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21;n+=1}') > t2; mv t2 12sol_ZY56_INS.merge.vcf
cat <(head -1 12sol_ZY56_INV.merge.vcf) <(grep -v '1CHROM' 12sol_ZY56_INV.merge.vcf |awk 'BEGIN{n=1;OFS="\t"}{print $1,$2,"SV_INV"n,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21;n+=1}') > t3; mv t3 12sol_ZY56_INV.merge.vcf
cat <(head -1 12sol_ZY56_TRANS.merge.vcf) <(grep -v '1CHROM' 12sol_ZY56_TRANS.merge.vcf |awk 'BEGIN{n=1;OFS="\t"}{print $1,$2,"SV_TRANS"n,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21;n+=1}') > t3; mv t3 12sol_ZY56_TRANS.merge.vcf
