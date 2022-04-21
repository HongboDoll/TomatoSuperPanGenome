#!/bin/bash

ls *fa | while read i
do
spe=`echo $i | awk -F '.' '{print $1}'`
./output_chr_fasta.py $i > ${spe}_chr.fa
done


for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65
do
qsub -N ${i} -l h_vmem=50g -q all.q -S /bin/bash -cwd -j y -V -b y "minimap2 -t 20 -ax asm5 --eqx ZY56_chr.fa ${i}_chr.fa | samtools view -b > ${i}_ZY56.bam"
samtools index ${i}_ZY56.bam
qsub -N ${i} -l h_vmem=50g -q all.q -S /bin/bash -cwd -j y -V -b y "/datalus/heqiang/biosoft/Python-3.5.5/bin/python3 /datalus/heqiang/biosoft/syri-master/syri/bin/syri -c ${i}_ZY56.bam -r ZY56_chr.fa -q ${i}_chr.fa -k -F B --prefix ${i}_ZY56 --nc 12"
done

for i in ZY59 Spenn ZY63 ZY64
do
qsub -N ${i} -l h_vmem=50g -q all.q -S /bin/bash -cwd -j y -V -b y "minimap2 -t 20  -ax asm10 --eqx ZY56_chr.fa ${i}_chr.fa | samtools view -b > ${i}_ZY56.bam"
#samtools index ${i}_ZY56.bam
qsub -N ${i} -l h_vmem=50g -q all.q -S /bin/bash -cwd -j y -V -b y "/datalus/heqiang/biosoft/Python-3.5.5/bin/python3 /datalus/heqiang/biosoft/syri-master/syri/bin/syri -c ${i}_ZY56.bam -r ZY56_chr.fa -q ${i}_chr.fa -k -F B --prefix ${i}_ZY56 --nc 12"
done


## INV
for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65 ZY59 Spenn ZY63 ZY64
do
cat ${i}_ZY56invOut.txt | grep '#' > ${i}_ZY56_INV.xls
        rm ${i}_ZY56_syri_INV.vcf
        cat ${i}_ZY56_INV.xls | while read n
        do
                chr=`echo $n | awk '{print $2}'`
                pos=`echo $n | awk '{print $3}'`
                ref_pos=`echo $n | awk '{print $2":"$3"-"$4}'`
                alt_pos=`echo $n | awk '{print $6":"$7"-"$8}'`
                ref_seq=`samtools faidx ZY56.fa $ref_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                alt_seq=`samtools faidx ${i}_chr.fa $alt_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                echo -e "${chr}\t${pos}\t.\t${ref_seq}\t${alt_seq}\t30\tPASS\t.\tGT\t1/1" >> ${i}_ZY56_syri_INV.vcf
        done
done

### DEL INS
for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65 ZY59 Spenn ZY63 ZY64
do
cat ${i}_ZY56syri.vcf|egrep 'INS|DEL'|awk 'length($4)-length($5)>50' > ${i}_ZY56_syri_DEL.vcf
cat ${i}_ZY56syri.vcf|egrep 'INS|DEL'|awk 'length($5)-length($4)>50' > ${i}_ZY56_syri_INS.vcf
done


## TRANS
for i in Heinz M82 ZY57 ZY58 ZY60 ZY61 ZY62 ZY65 ZY59 Spenn ZY63 ZY64
do
./output_translocation_from_vcf.py ${i}_ZY56syri.vcf > ${i}_ZY56_TRANS.xls
        rm ${i}_ZY56_syri_TRANS.vcf
        cat ${i}_ZY56_TRANS.xls | while read n
        do
                chr=`echo $n | awk '{print $2}'`
                pos=`echo $n | awk '{print $3}'`
                ref_pos=`echo $n | awk '{print $2":"$3"-"$4}'`
                alt_pos=`echo $n | awk '{print $6":"$7"-"$8}'`
                ref_seq=`samtools faidx ZY56.fa $ref_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                alt_seq=`samtools faidx ${i}_chr.fa $alt_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                echo -e "${chr}\t${pos}\t.\t${ref_seq}\t${alt_seq}\t30\tPASS\t.\tGT\t1/1" >> ${i}_ZY56_syri_TRANS.vcf
        done
done

