#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # merged.ont.v1.0.convert_GAL_filter300_${sv}.vcf


for line in i1:
	if '#' not in line:
		line = line.strip().split()
		if ',' not in line[4]:
			print('\t'.join(line[0:5]), '30', '.', '.', 'GT', sep='\t', end='\t')
		else:
			print('\t'.join(line[0:4]), line[4].split(',')[0], '30', '.', '.', 'GT', sep='\t', end='\t')
		for n in range(9, len(line)):
			k = line[n].split(':')[0]
			if n != (len(line) - 1):
				if k != '0/0' and k != './.' and k != '0/1' and k != '1/0':
					print('1/1', end='\t')
				else:
					print(k, end='\t')
			else:
				if k != '0/0' and k != './.' and k != '0/1' and k != '1/0':
					print('1/1')
				else:
					print(k)
	
