#!/bin/bash

ref=S_galapagense_canu_pilon.scaffolds.fasta

/vol1/agis/huangsanwen_group/lihongbo/bin/python /vol1/agis/huangsanwen_group/lihongbo/software/hic_pro/HiC-Pro_2.11.0-beta/bin/utils/digest_genome.py -r hindiii -o $ref.hindiii $ref

/vol1/agis/huangsanwen_group/lihongbo/bin/bowtie2-build ./$ref $ref

/vol1/agis/huangsanwen_group/lihongbo/software/faSize -detailed $ref > $ref.sizes

/vol1/agis/huangsanwen_group/lihongbo/software/hic_pro/HiC-Pro_2.11.0-beta/bin/HiC-Pro -i ./data -o hic-pro-result/ -c config-hicpro.txt

perl /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/PreprocessSAMs.pl /vol1/agis/huangsanwen_group/lihongbo/work/hic_pro_lachesis/hic-pro-result/bowtie_results/bwt2/sample1/S_galapagense_canu_pilon.bwt2pairs.bam $ref

         ###generate thousands of  running directories

perl /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/lachesisGreedy.pl 7 AAGCTT $ref $PWD S_galapagense_canu_pilon.bwt2pairs.REduced.paired_only.bam

         ### generate commdline

perl /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/runLachesis.pl run1/*ini

    ### merge  reports of all runs

while read line; do echo "perl /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/lachesisReport.pl 1 "$line"/REPORT.txt >"$line"_summary.txt && cat "$line"_summary.txt >> all_summary.txt" >> lachesis_summary.sh; done< run.list;


perl /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/convertHiCProMatrix2ordering.v2.pl ./main_results/ 15 500000 run23141 rawdata_500000_abs.bed rawdata_500000_iced.matrix
Rscript /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/heatmap.MWAH.R run23141.heatmap.txt run23141.chrom.breaks run23141.heatmap.png
perl /vol1/agis/huangsanwen_group/lihongbo/software/LACHESIS-master/CreateScaffoldedFasta.pl $ref  Lachesis/run23141