#!/bin/bash

source activate /home/jiayuxin/anaconda3/envs/ortho

rm -rf pep/OrthoFinder

orthofinder -f pep -a 52 -t 52 -T iqtree -M msa -ot

./convert_orthogroup_2_pan_matrix.py pep/OrthoFinder/Results*/Orthogroups/Orthogroups.tsv > tomato_matrix_core_dispensable_genes.cluster