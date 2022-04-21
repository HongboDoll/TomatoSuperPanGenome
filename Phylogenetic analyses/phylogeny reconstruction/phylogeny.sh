#!/bin/bash

source activate /home/jiayuxin/anaconda3/envs/ortho

orthofinder -f pep -a 52 -t 52 -T iqtree -M msa  -os

source /home/huyong/software/anaconda3/bin/activate collinerity

rm -rf muscle_out_pep; mkdir muscle_out_pep
ls $PWD/pep/OrthoFinder/Results_Mar06/Single_Copy_Orthologue_Sequences | while read i
do
	muscle -in $PWD/pep/OrthoFinder/Results_Mar06/Single_Copy_Orthologue_Sequences/${i} -out muscle_out_pep/${i}.out
	fa2phy.py -i muscle_out_pep/${i}.out -o muscle_out_pep/${i}.phylip
	awk 'NR!=1' muscle_out_pep/${i}.phylip | awk '{print $2}' > t && mv t muscle_out_pep/${i}.phylip
done

n=`ls muscle_out_pep/*.phylip | sed ':a;N;s/\n/ /g;ta'`
paste -d '' $n > muscle_out_pep/GR.phylip.all
awk 'BEIGN{n=1}{print ">"n"\n"$1;n+=1}' muscle_out_pep/GR.phylip.all > muscle_out_pep/GR.phylip.all.fa
n1=`grep -v '>' muscle_out_pep/GR.phylip.all.fa | wc -l | awk '{print $1}'`
n2=`sed 's/ //g' muscle_out_pep/GR.phylip.all.fa | grep -v '>' | head -1 |wc | awk '{print $3-1}'`
echo -e "$n1 $n2" > muscle_out_pep/head
cat muscle_out_pep/head <(paste -d '\t' first_col  <(grep -v '>' muscle_out_pep/GR.phylip.all.fa)) | sed 's/\t/            /g' > single_copy_gene_family_pep_4_tree.phylip

phyml -i single_copy_gene_family_pep_4_tree.phylip -d aa -b 100 --model JTT -f e -v 0.576 -a 0.886 --nclasses 4 --search SPR -t e 

mv proteins.phy_phyml_tree.txt tomato_single_copy_phyml_tree.nwk
