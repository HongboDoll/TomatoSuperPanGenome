#!/bin/bash

##### reroot the tree using Figtree

# tomato_single_copy_phyml_tree.nwk

#### remove branch lengths and bootstrap values

java -jar $li/bin/PareTree1.0.2.jar -t O -f tomato_single_copy_phyml_tree_rooted.nwk -nbs -topo

num_spe=`cat tomato_single_copy_phyml_tree_rooted_nbs.nwk | awk -F ',' '{print NF}'`
cat <(echo -e "$num_spe 1") tomato_single_copy_phyml_tree_rooted_nbs.nwk > t && mv t tomato_single_copy_phyml_tree_rooted.nwk

tree=tomato_single_copy_phyml_tree_rooted.nwk

#### add fossil calibration time for the two tree files, as well as the constrait for the root node
cp $tree ${tree}_baseml
cp $tree ${tree}_mcmctree

#### baseml

source activate /home/jiayuxin/anaconda3/envs/paml

baseml baseml.ctl


#### mcmctree

num_spe=`cat tomato_single_copy_phyml_tree_rooted.nwk | awk -F ',' '{print NF}'`
cat <(echo -e "$num_spe 1") ${tree}_mcmctree > t && mv t ${tree}_mcmctree

#### Prepare the .ctl file for mcmctree manually
## seqfile = the alignment file used in last step (baseml)
## treefile = modify the tree, add fossil divergence time for calibration
## IMPORTANT: rgene_gamma = 1 / (subsititution_rate)
############################################### 
#
#### First run, use parameter 'usedata = 3' to generate 'out.BV'
mcmctree mcmctree.usedata3.ctl &>log.first
#
#### rename
cp out.BV in.BV
#
#### Second run, use parameter 'usedata = 2'
mcmctree mcmctree.usedata2.ctl &>log.second

