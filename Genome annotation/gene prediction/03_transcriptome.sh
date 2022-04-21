#!/bin/bash

ref=S_galapagense_canu_pilon.chr.fasta
spe=S_gal
threads=18

hisat2-build -p 20 $ref re1
hisat2 -x re1 --dta -p 20 -1 CUChhxTERAAPE_1.clean.fq.gz -2 CUChhxTERAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_stem.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTHRAAPE_1.clean.fq.gz -2 CUChhxTHRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_female.flower.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTBRAAPE_1.clean.fq.gz -2 CUChhxTBRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_Expanded.ovary.Fertilized.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTFRAAPE_1.clean.fq.gz -2 CUChhxTFRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_leaf.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTGRAAPE_1.clean.fq.gz -2 CUChhxTGRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_male.flower.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTIRAAPE_1.clean.fq.gz -2 CUChhxTIRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_tendril.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTDRAAPE_1.clean.fq.gz -2 CUChhxTDRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_root.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTCRAAPE_1.clean.fq.gz -2 CUChhxTCRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_Expanded.ovary.Unfertilized.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTARAAPE_1.clean.fq.gz -2 CUChhxTARAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_ovary.bam
hisat2 -x re1 --dta -p 20 -1 CUChhxTJRAAPE_1.clean.fq.gz -2 CUChhxTJRAAPE_2.clean.fq.gz | samtools view -bS -@ 20 - > ${spe}_base.of.tendril.bam

samtools merge -@ 20 ${spe}_10tissues.bam ${spe}_stem.bam ${spe}_female.flower.bam ${spe}_Expanded.ovary.Fertilized.bam ${spe}_leaf.bam ${spe}_male.flower.bam ${spe}_tendril.bam ${spe}_root.bam ${spe}_Expanded.ovary.Unfertilized.bam ${spe}_ovary.bam ${spe}_base.of.tendril.bam && rm ${spe}_stem.bam ${spe}_female.flower.bam ${spe}_Expanded.ovary.Fertilized.bam ${spe}_leaf.bam ${spe}_male.flower.bam ${spe}_tendril.bam ${spe}_root.bam ${spe}_Expanded.ovary.Unfertilized.bam ${spe}_ovary.bam ${spe}_base.of.tendril.bam
samtools sort -@ 20 -m 2G ${spe}_10tissues.bam -o ${spe}_10tissues.sort.bam && rm  ${spe}_10tissues.bam

stringtie -p 20 ${spe}_10tissues.sort.bam -o stringtie_out

sed 's/StringTie/Cufflinks/g' stringtie_out > transcripts.gtf

/share/fg2/lihb/software/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl stringtie_out $ref > transcripts.fasta
/share/fg2/lihb/software/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl stringtie_out > transcripts.gff3
/share/fg2/lihb/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcripts.fasta -m 100

#makeblastdb -in uniprot_sprot_plants.fa -dbtype prot
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db uniprot_sprot_plants.fa -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 20 > blastp.outfmt6

hmmsearch --cpu 18 -o ttt --domtblout hmmsearch.tmp Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep
awk 'BEGIN{OFS=FS=" "} NR<=3{print}; NR>3{tmp=$1; $1=$4; $4=tmp; tmp=$2; $2=$5; $5=tmp; print}' hmmsearch.tmp > pfam.domtblout

/share/fg2/lihb/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t transcripts.fasta --single_best_only --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

/share/fg2/lihb/software/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

