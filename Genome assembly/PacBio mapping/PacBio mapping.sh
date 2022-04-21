#!/bin/bash

for i in S_galapagense S_lycopersicum_var_cerasiforme S_lycopersicum_M82 #S_chilense S_corneliomulleri S_habrochaites S_lycopersicoides S_neorickii S_peruvianum
do
       qsub -N ${i}_pacbio -q all.q -l h_vmem=50g -S /bin/bash -cwd -j y -V -b y "minimap2 -t 10 -ax map-pb 13_tomato_genomes/${i}.fa Canu_corected/${i}*.fasta.gz | samtools sort -@ 10 > ${i}.pacbio.sort.bam"
       qsub -N ${i}_index -q all.q -l h_vmem=10g  -S /bin/bash -cwd -j y -V -b y "samtools index ${i}.pacbio.sort.bam"
size=`seq_n50.pl 13_tomato_genomes/${i}.fa | grep Total | awk '{print $2}'`
for n in {01..12}
do
	qsub -N ${i}_depth -q all.q -l h_vmem=20g -S /bin/bash -cwd -j y -V -b y "bamCoverage -b ${i}.pacbio.sort.bam -o ${i}.Sly${n}.pacbio.bedgraph -of bedgraph --binSize 200000 --normalizeUsing RPGC --effectiveGenomeSize $size --smoothLength 300000 --region Sly${n} -p 10"
done
cat ${i}.Sly*pacbio.bedgraph > ${i}.pacbio.bedgraph
./plot_depth.R ${i}.pacbio.bedgraph ${i}.pacbio.bedgraph.pdf

done
