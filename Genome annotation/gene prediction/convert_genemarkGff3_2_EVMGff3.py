#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # ${i}genemark.gff3

for line in i1:
	line = line.strip().split()
	if line[2] == 'exon':
		print('\t'.join(line[0:8]), end='\t')
		Id = line[8].split(';')[0].split('=')[1]
		parent = line[8].split(';')[1].split('=')[1]
		print('ID='+parent+'.'+line[8].split(';')[0].split('=')[1]+';'+line[8].split(';')[1]+';')
	elif line[2] == 'CDS':
		print('\t'.join(line[0:8]), end='\t')
		Id = line[8].split(';')[0].split('=')[1]
		parent = line[8].split(';')[1].split('=')[1]
		print(line[8].split(';')[0]+'.'+parent+';'+line[8].split(';')[1]+';')
	else:
		print('\t'.join(line[0:]))
		

