#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 12sol_ZY56_SV_syri_merge_INS.vcf
i2 = open(sys.argv[2])  # 12sol_ZY56_SV_syri_merge_INV_no_wild.vcf

sv = {}
for line in i2:
    if '#' not in line:
        line = line.strip().split()
        sv[line[2]] = ''
        
for line in i1:
    if '#' in line:
        print(line.strip())
    else:
        line = line.strip().split()
        if line[2] not in sv:
            print('\t'.join(line))

