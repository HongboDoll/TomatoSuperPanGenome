#!/bin/bash

######### extract four classes of SVs
for sv in CNV INV  INS DEL
do
	grep $sv 12sol_ZY56_SV_syri_merge.vcf > 12sol_ZY56_SV_syri_merge_${sv}.vcf
done

for sv in CNV INV  INS DEL
do
	if [ $sv = "CNV" ];then
		n='DUP'
		grep "SVTYPE=$n" merged.ont.v1.0.convert_GAL_filter300.vcf > merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf
		grep "SVTYPE=$n" merged.ont.v1.0.vcf > merged.ont.v1.0_${sv}.vcf
	else
		grep "SVTYPE=$sv" merged.ont.v1.0.convert_GAL_filter300.vcf > merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf
		grep "SVTYPE=$sv" merged.ont.v1.0.vcf > merged.ont.v1.0_${sv}.vcf
	fi
done


######## revise INV and DUP in cell100 (SURVIVOR does not support implicit ALT)
awk '$5!~/>/' merged.ont.v1.0.convert_GAL_filter300_INS.vcf > t && mv t merged.ont.v1.0.convert_GAL_filter300_INS.vcf
awk '$5!~/>/' merged.ont.v1.0.convert_GAL_filter300_DEL.vcf > t && mv t merged.ont.v1.0.convert_GAL_filter300_DEL.vcf
awk '$5=="<INV>"' merged.ont.v1.0.convert_GAL_filter300_INV.vcf | sed 's/<INV>/CTCTCTCTCTCTCTCT/g' > t && mv t merged.ont.v1.0.convert_GAL_filter300_INV.vcf
awk '$5=="<INV>"' merged.ont.v1.0_INV.vcf | sed 's/<INV>/CTCTCTCTCTCTCTCT/g' > t && mv t merged.ont.v1.0_INV.vcf
awk '$5=="<DUP>"' merged.ont.v1.0.convert_GAL_filter300_CNV.vcf | sed 's/<DUP>/CTCTCTCTCTCTCTCT/g' > t && mv t merged.ont.v1.0.convert_GAL_filter300_CNV.vcf
awk '$5=="<DUP>"' merged.ont.v1.0_CNV.vcf | sed 's/<DUP>/CTCTCTCTCTCTCTCT/g' > t && mv t merged.ont.v1.0_CNV.vcf


######## convert format, remove INFO column
for sv in CNV INV  INS DEL
do
   cat 100_cell_head <(./convert_12_sol_100_cell_vcf.py merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf) > t; mv t merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf
   cat 100_cell_head <(./convert_12_sol_100_cell_vcf.py merged.ont.v1.0_${sv}.vcf) > t; mv t merged.ont.v1.0_${sv}.vcf
   cat 12_sol_head <(./convert_12_sol_100_cell_vcf.py 12sol_ZY56_SV_syri_merge_${sv}.vcf) > t; mv t 12sol_ZY56_SV_syri_merge_${sv}.vcf
done

#

######### handle those close variants (shorter than 50 bp flanking)
for sv in CNV INV  INS DEL
do
	ls merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf > vcf_4_merge
	SURVIVOR merge vcf_4_merge 50 1 0 0 0 0 merged.vcf
	./revise_vcf_GT1.py merged.vcf merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf > t && mv t merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf
	./revise_merge_vcf_CHROM.py merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf head_100 > t && mv t merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf
	ls 12sol_ZY56_SV_syri_merge_${sv}.vcf > vcf_4_merge
	SURVIVOR merge vcf_4_merge 50 1 0 0 0 0 merged.vcf
	./revise_vcf_GT1.py merged.vcf 12sol_ZY56_SV_syri_merge_${sv}.vcf > t && mv t 12sol_ZY56_SV_syri_merge_${sv}.vcf 
	./revise_merge_vcf_CHROM.py 12sol_ZY56_SV_syri_merge_${sv}.vcf head_12 > t && mv t 12sol_ZY56_SV_syri_merge_${sv}.vcf
done

######## merge the two sets

grep 'CHROM' 12sol_ZY56_SV_syri_merge.vcf|cut -f '10-999' > head_12
grep 'CHROM' merged.ont.v1.0.convert_GAL_filter300.vcf |cut -f '10-999' > head_100

for sv in CNV INV  INS DEL
do
	ls 12sol_ZY56_SV_syri_merge_${sv}.vcf merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf > vcf_4_merge
	SURVIVOR merge vcf_4_merge 50 1 0 0 0 0 merged.vcf
	./revise_vcf_GT2.py merged.vcf 12sol_ZY56_SV_syri_merge_${sv}.vcf merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf head_12 head_100 > 12_100_merged_ZY56_${sv}.vcf
done


######## check VCF REF sequences Versus reference

rm -rf revise_ref && mkdir revise_ref
for sv in INS DEL INV CNV
do
	for i in `seq 1 400`
	do
       ./output_work2_shell_${sv}.py $i > work2_${sv}_${i}.sh && chmod 755 work2_${sv}_${i}.sh
       #sbatch -J ${sv}_${i} -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1  -e %x.err -o %x.out "./work2_${sv}_${i}.sh"
	done
done

for sv in INS INV
do
    for i in `seq 1 400`
    do
		sbatch -J ${sv}_${i} -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1  -e %x.err -o %x.out "./work2_${sv}_${i}.sh"
	done
done

for sv in DEL CNV
do
    for i in `seq 1 400`
    do
        sbatch -J ${sv}_${i} -p low -N 1 --ntasks-per-node=1  -e %x.err -o %x.out "./work2_${sv}_${i}.sh"
    done
done
                                                                                                                                                                          
for sv in INS DEL INV CNV
do
echo -e """#!/bin/bash  
 cat <(grep '#' 12_100_merged_ZY56_INS.vcf) <(cat revise_ref/*_${sv} | sort -k1,1 -k2,2n) > 12_100_merged_ZY56_${sv}_revised.vcf 
""" > ${sv}_merge.sh && chmod 755 ${sv}_merge.sh
done

for sv in INS DEL INV CNV
do
	sbatch -J ${sv}_merge -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1 -e %x.err -o %x.out "./${sv}_merge.sh"
done

### sv in wild species

for sv in INS DEL INV CNV
do
	./filter_no_wild_SV.py 12sol_ZY56_SV_syri_merge_${sv}.vcf 12sol_ZY56_SV_syri_merge_${sv}_no_wild.vcf > 12sol_ZY56_SV_syri_merge_${sv}_wild.vcf
done

for sv in INS DEL INV
do
    ls 12sol_ZY56_SV_syri_merge_${sv}_wild.vcf merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf > vcf_4_merge
    SURVIVOR merge vcf_4_merge 50 1 0 0 0 0 merged_${sv}_wild.vcf
	grep -v '#' 12sol_ZY56_SV_syri_merge_${sv}_wild.vcf | wc -l
	grep 'SUPP=2' merged_${sv}_wild.vcf | wc -l
	
done

for i in DEL INS INV; do grep -v '#' merged_${i}_wild.vcf|grep -v 'SUPP=2' > merged_${i}_wild_our_unique.vcf; done

rm merged_INS_DEL_wild_our_unique_impacted_genes.xls
for sv in DEL INS
do
	java -jar /public/agis/huangsanwen_group/lihongbo/software/snpEff/snpEff.jar ZY56 -ud 2000 merged_${sv}_wild_our_unique.vcf > merged_${sv}_wild_our_unique.vcf.snpeff
	./output_exon_gene.py merged_${sv}_wild_our_unique.vcf.snpeff | sort | uniq >> merged_INS_DEL_wild_our_unique_impacted_genes.xls
done

