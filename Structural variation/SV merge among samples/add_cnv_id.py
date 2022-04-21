#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # tmp 
i2 = open(sys.argv[2])  # 12sol_ZY56_CNV.merge.vcf

d = {}
for line in i1:
	line = line.strip().split()
	if 'Q' in line[3]:
		d[line[0]+'~'+line[1]] = 'CNV-Q'
	else:
		d[line[0]+'~'+line[1]] = 'CNV-R'

n = 0
for line in i2:
	if '1CHROM' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		if line[0]+'~'+line[1] in d:
			print(line[0], line[1], 'SV_'+d[line[0]+'~'+line[1]]+str(n), '\t'.join(line[3:]), sep='\t')
		else:
			print(line[0], line[1], 'SV_'+'CNV'+str(n), '\t'.join(line[3:]), sep='\t')
	n += 1
