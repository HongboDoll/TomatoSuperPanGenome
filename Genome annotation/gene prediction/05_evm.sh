#!/bin/bash

##### EVM main

spe=S_gal
ref=S_galapagense_canu_pilon.chr.fasta

cat ${spe}.snap.gff3 ${spe}.augustus.gff3 ${spe}.genemark.gff3 <(sed 's/transdecoder/OTHER_PREDICTION/g' transcripts.fasta.transdecoder.genome.gff3) > evm_denovo.gff3
cat test_${spe}.pasa_assemblies.gff3 > evm_rnaseq.gff3

EVM_HOME='/share/fg2/lihb/software/EVidenceModeler-1.1.1/'

rm -rf tig*
$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome $ref --gene_predictions evm_denovo.gff3 --protein_alignments ${spe}.exonerate.gff3 --transcript_alignments evm_rnaseq.gff3 --segmentSize 800000 --overlapSize 20000 --partition_listing partitions_list.out
$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome $ref --weights $PWD/weights.txt --gene_predictions evm_denovo.gff3 --protein_alignments ${spe}.exonerate.gff3 --transcript_alignments evm_rnaseq.gff3 --output_file_name evm.out --partitions partitions_list.out > commands.list

split -a 2 -d -l 70 commands.list commands_
for i in {00..20}
do
mv commands_$i commands_${i}.sh; chmod 755 commands_${i}.sh
nohup ./commands_${i}.sh &
done
###########################
#
EVM_HOME='/share/fg2/lihb/software/EVidenceModeler-1.1.1/'
$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome $ref
cat */evm.out.gff3 > ${spe}.evm_out.gff3

