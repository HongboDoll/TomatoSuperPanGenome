#!/bin/bash

ref=S_galapagense_canu_pilon.contigs.fasta

/vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/RefAligner/1.0/RefAligner -i S_gal_RawMolecules.bnx -minlen 100 -merge -o S_gal_RawMolecules.filter100kb -bnx

/usr/bin/perl /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/HybridScaffold/12162019/scripts/fa2cmap_multi_color.pl -i $ref -e CTTAAG 1

python /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/Pipeline/12162019/align_bnx_to_cmap.py --prefix AM --mol S_gal_RawMolecules.filter100kb.bnx --ra /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/RefAligner/1.0/ --ref assembly_CTTAAG_0kb_0labels.cmap --nthreads 20 --output S_gal_bnx2cmap --snrFilter 2 --color 1

python /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/Pipeline/12162019/pipelineCL.py -T 20 -N 4 -f 0 -i 5 -y -b S_gal_RawMolecules.filter100kb.bnx -l S_gal_BN_asm.output -t /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/RefAligner/1.0/ -a /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/RefAligner/1.0/optArguments_nonhaplotype_noES_noCut_BG_saphyr.xml -R

/usr/bin/perl /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/HybridScaffold/1.0/hybridScaffold.pl -n assembly.fasta -b kbs-mac-74_bng_contigs2017.cmap -c /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/HybridScaffold/1.0/hybridScaffold_config.xml -r /vol1/agis/huangsanwen_group/lihongbo/software/Solve3.5_12162019/RefAligner/1.0/RefAligner -o S_gal_hybrid -B 2 -N 2 -f    
