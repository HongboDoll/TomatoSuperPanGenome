#!/bin/bash                                                                                                                                      

ref=S_galapagense_canu_pilon.chr.fasta
spe=S_gal
thread=10

############ trinity without ref  #### 
#
/share/fg2/lihb/software/trinityrnaseq-v2.10.0/Trinity --seqType fq --max_memory 100G --left s80_combine_1.clean.fq.gz,9930_39_tissues_1.fq.gz --right s80_combine_2.clean.fq.gz,9930_39_tissues_2.fq.gz -output ./S_gal_trinity_out_no_ref --min_kmer_cov 2 --trimmomatic --normalize_reads --CPU $thread 
#
############### trinity with ref

/share/fg2/lihb/software/trinityrnaseq-v2.10.0/Trinity --genome_guided_bam ./${spe}_10tissues.sort.bam  --max_memory 50G --genome_guided_max_intron 10000 --output ./S_gal_trinity_out_with_ref --CPU $thread

################ PASA

cat ./S_gal_trinity_out_no_ref/Trinity.fasta ./S_gal_trinity_out_with_ref/Trinity-GG.fasta > transcripts.fasta
cat transcripts.denovo.fasta ./S_gal_trinity_out_with_ref/Trinity-GG.fasta > transcripts.fasta
/share/fg2/lihb/software/PASApipeline.v2.4.1/misc_utilities/accession_extractor.pl < transcripts.fasta > tdn.accs
/share/fg2/lihb/software/PASApipeline.v2.4.1/seqclean/seqclean/seqclean  transcripts.fasta
/share/fg2/lihb/software/PASApipeline.v2.4.1/scripts/Launch_PASA_pipeline.pl -c /share/fg2/lihb/software/PASApipeline.v2.4.1/pasa_conf/pasa.alignAssembly.Template.txt --trans_gtf transcripts.gtf --TDN tdn.accs -C  -R -g $ref -t transcripts.fasta.clean -T -u transcripts.fasta --ALIGNERS blat --CPU $thread
/share/fg2/lihb/software/PASApipeline.v2.4.1/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta test_80.assemblies.fasta --pasa_transcripts_gff3 test_80.pasa_assemblies.gff3
/share/fg2/lihb/software/PASApipeline.v2.4.1/scripts/pasa_asmbls_to_training_set.extract_reference_orfs.pl  test_80.assemblies.fasta.transdecoder.genome.gff3 > best_candidates.gff3

