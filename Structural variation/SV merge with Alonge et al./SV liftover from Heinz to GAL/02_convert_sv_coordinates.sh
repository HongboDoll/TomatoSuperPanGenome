#!/bin/bash

rm -rf convert_vcf1 && mkdir convert_vcf1
rm -rf tmp && mkdir tmp
for i in `seq 1 400`
do
	./output_work2_shell.py $i > work2_${i}.sh && chmod 755 work2_${i}.sh
	sbatch -J 2_${i} -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1  -e %x.err -o %x.out "./work2_${i}.sh"
done

cat convert_vcf1/*vcf  | ./merge_vcf.py merged.ont.v1.0.vcf > merged.ont.v1.0.convert_GAL_filter300.vcf

cat <(grep '#' merged.ont.v1.0.convert_GAL_filter300.vcf) <(grep -v '#' merged.ont.v1.0.convert_GAL_filter300.vcf | sort -k1,1 -k2,2n) > t

mv t merged.ont.v1.0.convert_GAL_filter300.vcf

