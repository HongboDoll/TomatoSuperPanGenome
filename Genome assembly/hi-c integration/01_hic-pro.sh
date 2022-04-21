#!/bin/bash

for spe in S_chilense S_chmielewskii S_galapagense S_lycopersicum_M82 S_lycopersicum_var_cerasiforme S_pimpinellifolium S_corneliomulleri S_habrochaites S_lycopersicoides S_neorickii S_peruvianum
do

rm -rf $spe; mkdir -p ${spe}/sample; cd ${spe}/sample
ln -s $PWD/../../hi-c_data/${spe}/*gz .
cd -

rm -rf /datalus/heqiang/tomato/2022_1_15_revision/heterozygosity/hi-c_mapping/${spe}_ref; mkdir /datalus/heqiang/tomato/2022_1_15_revision/heterozygosity/hi-c_mapping/${spe}_ref; 
cd /datalus/heqiang/tomato/2022_1_15_revision/heterozygosity/hi-c_mapping/${spe}_ref

../retain_only_chr.py ../13_tomato_genomes/${spe}.fa > t
seqkit sort -n t > ../13_tomato_genomes/${spe}_chr_sort.fa && rm t

ln -s ../13_tomato_genomes/${spe}_chr_sort.fa ${spe}.fa
cd -
echo -e """#!/bin/bash
cd ${spe}_ref
bowtie2-build ${spe}.fa ${spe}.fa

out_len.py ${spe}.fa | sed 's/>//g' > ${spe}.fa.sizes

/data2/liujun/HQ/programs/anaconda2/bin/python /data2/liujun/HQ/programs/Hic_pro/HiC-Pro_2.11.1/bin/utils/digest_genome.py -r hindiii -o ${spe}.fa.hindiii ${spe}.fa
cd -

sed "s/ZY56/${spe}/g" config.txt > ${spe}_config.txt
rm -rf  ./${spe}_pro
source activate HiC-Pro_v3.1.0
/datalus/heqiang/biosoft/hicpro/HiC-Pro_3.1.0/bin/HiC-Pro -i ./${spe} -o ./${spe}_pro -c ${spe}_config.txt
""" > ${spe}_run.sh && chmod 755 ${spe}_run.sh
done


for spe in S_chilense S_chmielewskii S_galapagense S_lycopersicum_M82 S_lycopersicum_var_cerasiforme S_pimpinellifolium S_corneliomulleri S_habrochaites S_lycopersicoides S_neorickii S_peruvianum
do
qsub -N ${spe}_hic -q all.q -l h_vmem=50g -S /bin/bash -cwd -j y -V -b y "./${spe}_run.sh"
done

win_size=100000
#
for spe in S_chilense S_chmielewskii S_galapagense S_lycopersicum_M82 S_lycopersicum_var_cerasiforme S_pimpinellifolium S_corneliomulleri S_habrochaites S_lycopersicoides S_neorickii S_peruvianum
do
#### whole genome heat maps
/data2/liujun/HQ/programs/anaconda2/bin/python /data2/liujun/HQ/programs/HiCPlotter/HiCPlotter.py -tri 1 -f ${spe}_pro/hic_results/matrix/sample/iced/${win_size}/sample_${win_size}_iced.matrix -bed ${spe}_pro/hic_results/matrix/sample/raw/${win_size}/sample_${win_size}_abs.bed -wg 1 -n ${spe}_whole_genome -fh 0 -r $win_size -hmc 3 -o ${spe}_HiCPlot_${win_size}_whole_genome -ext pdf -dpi 300 # -mm 15

### intra-chromosomal heat maps
for i in {01..12}
do
	/data2/liujun/HQ/programs/anaconda2/bin/python /data2/liujun/HQ/programs/HiCPlotter/HiCPlotter.py -tri 1 -f ${spe}_pro/hic_results/matrix/sample/iced/${win_size}/sample_${win_size}_iced.matrix -bed ${spe}_pro/hic_results/matrix/sample/raw/${win_size}/sample_${win_size}_abs.bed -wg 0 -chr Sly${i} -n Sly${i}  -fh 0 -r $win_size -hmc 3 -o ${spe}_HiCPlot_${win_size}_Sly${i}  -ext pdf -dpi 300 # -mm 15
done
done


