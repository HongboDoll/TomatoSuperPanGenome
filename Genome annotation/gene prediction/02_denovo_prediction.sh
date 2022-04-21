#!/bin/bash

ref=S_galapagense_canu_pilon.chr.fasta


############### Genescan

genscan  genscan_smat  $ref > ${ref}.predict

############# augustus set train
##

/share/fg2/lihb/software/augustus-3.1/scripts/gff2gbSmallDNA.pl best_candidates.gff3 $ref 1000 genes.raw.gb
/share/fg2/lihb/software/augustus-3.1/scripts/new_species.pl --species=S_gal_for_bad_removingall${spe}
/share/fg2/lihb/software/augustus-3.1/bin/etraining --species=S_gal_for_bad_removingall${spe} --/genbank/verbosity=2 genes.raw.gb  > train.err
cat train.err | perl -pe 's/.*in sequence (\S+):/$1/' > badgenes.lst
/share/fg2/lihb/software/augustus-3.1/scripts/filterGenes.pl  badgenes.lst genes.raw.gb >train.gb
/share/fg2/lihb/software/augustus-3.1/scripts/randomSplit.pl train.gb 100
/share/fg2/lihb/software/augustus-3.1/scripts/new_species.pl --species=S_galall${spe}
/share/fg2/lihb/software/augustus-3.1/bin/etraining --species=S_galall${spe} train.gb.train > train.out
/share/fg2/lihb/software/augustus-3.1/bin/augustus --species=S_galall${spe} train.gb.test |tee firsttest.out
/share/fg2/lihb/software/augustus-3.1/bin/etraining --species=S_galall${spe} train.gb.train
/share/fg2/lihb/software/augustus-3.1/bin/augustus --species=S_galall${spe} train.gb.test |tee second.out

############ run augustus
ref=S_galapagense_canu_pilon.chr.fasta

seqkit split -p 20 $ref
rm -rf aug_out; mkdir aug_out
for i in {01..20}
do
nohup /share/fg3/lihb/software/augustus-3.1/bin/augustus --AUGUSTUS_CONFIG_PATH=/share/fg3/lihb/software/augustus-3.1/config --species=S_galall${spe} ${ref}.split/*part_0${i}* > aug_out/${i}.out &
done

cat aug_out/*out > ${ref}.out

grep -v "^#" ${ref}.out > gy14.augustus.format.out
/share/fg2/lihb/software/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl ${spe}.augustus.format.out >  ${spe}.augustus.gff3

######### glimmerhmm
awk '$3=="OK" {print$2}' validate.log > ok_name

python /vol1/agis/huangsanwen_group/lihongbo/work/pan_genome/80_S_gal_annotaion/scripts/extrac_ok_gene.py ok_name best_candidates_format.gff3 best_candidates_format_ok
python /vol1/agis/huangsanwen_group/lihongbo/work/pan_genome/80_S_gal_annotaion/scripts/eachgff2zff.py best_candidates_format_ok > best_candidates_format_ok.zff
python /vol1/agis/huangsanwen_group/lihongbo/work/pan_genome/80_S_gal_annotaion/scripts/format_trans.py best_candidates_format_ok.zff exon_ok.fa


awk '{if($1~/^>..../)print}' best_candidates_format_ok.zff >1
uniq 1 >pos
rm 1
python /vol1/agis/huangsanwen_group/lihongbo/work/pan_genome/80_S_gal_annotaion/scripts/change_pos_snap.py genome.dna pos 1
grep -v "^$" 1 >genome_ok.dna
rm 1
/public/agis/huangsanwen_group/liqing/software/GlimmerHMM/train/trainGlimmerHMM genome_ok.dna  exon_ok.fa -d S_gal_ok

############# glimmerhmm
#
grep '>' $ref > pos
python  /vol1/agis/huangsanwen_group/lihongbo/work/pan_genome/80_S_gal_annotaion/scripts/extract_each_seq.py pos $ref extract_each_seq_command
split -a 1 -d -l 200 extract_each_seq_command extract_each_seq_command_
mkdir seq
for i in {0..12}
do
               mv extract_each_seq_command_$i extract_each_seq_command_$i.sh
               chmod 755 extract_each_seq_command_$i.sh
               quick_qsub =extract= {-q queue2} ./extract_each_seq_command_$i.sh
done

########################################################################################

./split_glimmer.py $n_contig $spe > glimmerhmm_  # $n_contig is the number of contig of ref
split -a 2 -d -l 100 glimmerhmm_ glimmerhmm.sh_
mkdir glimmerhmm
for i in {00..10}
       do
       chmod 755 glimmerhmm.sh_${i}
       quick_qsub =${i}_glimmer= {-q queue1 -l nodes=1:ppn=1} ./glimmerhmm.sh_${i}
       done

#########################################################################################

cat ./glimmerhmm/*.out > glimmerhmm.out

grep -v "^$" glimmerhmm.out > glimerhmm_out_format
./glimmer_format.py glimerhmm_out_format glimerhmm_evm.gff3
