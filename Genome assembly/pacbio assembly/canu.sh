#!/bin/bash

canu -p S_galapagense -d S_galapagense_canu genomeSize=900m -pacbio-raw S_galapagense.pacbio.fastq.gz
awk '{if($1~/>/){print $1}else{print}}' S_galapagense_canu.contigs.fasta > tmp && mv tmp S_galapagense_canu.contigs.fasta



