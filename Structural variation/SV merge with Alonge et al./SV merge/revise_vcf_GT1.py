#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # merged.vcf without genotype
i2 = open(sys.argv[2])  # original.vcf with genotype

geno = {}
for line in i2:
	if '#' not in line:
		line = line.strip().split()
		if line[0]+'_'+line[1] not in geno:
			geno[line[0]+'_'+line[1]] = line[9:]
		else:
			for n in range(0, len(geno[line[0]+'_'+line[1]])):
				if line[n+9] != '0/0' and line[n+9] != './.' and (geno[line[0]+'_'+line[1]][n] == '0/0' or geno[line[0]+'_'+line[1]][n] == './.'):
					geno[line[0]+'_'+line[1]][n] = line[n+9]
		for k in geno[line[0]+'_'+line[1]]:
			if k != '0/0' and k != './.' and k != '0/1' and k != '1/0':
				k = '1/1'
				

for line in i1:
	if '#' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		sv = line[-1].split(':')[-1].split(',')
		gt = []
		if len(sv) > 1:
			for k in sv:
				if k.split('-')[0] in geno:
					if not gt:
						gt = geno[k.split('-')[0]]
					else:
						for n in range(0, len(gt)):
							if geno[k.split('-')[0]][n] != '0/0' and geno[k.split('-')[0]][n] != './.' and (gt[n] == '0/0' or gt[n] == './.'):
								gt[n] = geno[k.split('-')[0]][n]
			for k in gt:
				if k != '0/0' and k != './.' and k != '0/1' and k != '1/0':
					k = '1/1'
			print('\t'.join(line[0:6]), '.', '.', 'GT', '\t'.join(gt), sep='\t')
		else:
			print('\t'.join(line[0:6]), '.', '.', 'GT', '\t'.join(geno[line[0]+'_'+line[1]]), sep='\t')

