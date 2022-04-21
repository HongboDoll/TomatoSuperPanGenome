#!/bin/bash

./output_go_ipr_from_annot.py S_galapagense.pos.go.ipr.txt S_galapagense_GO.xls2 S_galapagense_IPR.xls > S_galapagense_GO.xls

./go.R S_galapagense_GO.xls2 merged_INS_DEL_wild_our_unique_impacted_genes.xls wild_SV_GO_ratio_infor.xls

