#!/bin/bash

spe=S_gal
ref=S_galapagense_canu_pilon.chr.fasta
threads=52
out=${spe}_maker
rm -rf tmp; mkdir tmp
rm -rf ${spe}_maker.maker.output

for i in Arabidopsis_thaliana Oryza_sativa Solanum_lycopersicum Solanum_tuberosum
do
	gff3=${i}_gene.gff3
	ref1=${i}.fa
	mkdir $PWD/${i}_outdir
	java -jar GeMoMa-1.8.jar CLI GeMoMaPipeline threads=$threads outdir=$PWD/${i}_outdir GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$ref i=$spe a=$gff3 g=$ref1
done



