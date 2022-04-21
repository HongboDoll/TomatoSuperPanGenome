#!/bin/bash

source /home/huyong/software/anaconda3/bin/activate /home/huyong/software/anaconda3/envs/busco-5.0

busco -i S_galapagense_canu.chr.fasta -l /public/agis/huangsanwen_group/lihongbo/database/embryophyta_odb10 -o S_galapagense_canu.chr.fasta_busco -m geno -c 52 --offline -f --augustus