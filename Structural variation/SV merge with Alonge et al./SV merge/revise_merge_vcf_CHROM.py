#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # merged.vcf without genotype
i2 = open(sys.argv[2])  # head_12

acc1 = []
for line in i2:
	acc1 = line.strip().split()

for line in i1:
	if '#' in line and '#CHROM' not in line:
		print(line.strip())
	elif '#' in line and '#CHROM' in line:
		line = line.strip().split()
		print('\t'.join(line[0:9]), end='\t')
		for n in range(0, len(acc1)):
			if n != (len(acc1) - 1):
				print(acc1[n], end='\t')
			else:
				print(acc1[n])
	else:
		print(line.strip())

