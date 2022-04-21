#!/usr/bin/env python3

import sys,re

i1 = open(sys.argv[1])  # vcf

for line in i1:
	if '#' not in line:
		line = line.strip().split()
		#r1 = re.search('\|upstream_gene_variant\|MODIFIER\|(\w+)\|(\w+)\|transcript\|', line[7]) # |upstream_gene_variant|MODIFIER|Sgal01g000130|
		r2 = re.search('\|MODERATE\|(\w+)\|(\w+)\|transcript\|', line[7])
		r3 = re.search('\|HIGH\|(\w+)\|(\w+)\|transcript\|', line[7])
		
		#if r1:
	#		print(r1.group(1))
		if r2:
			print(r2.group(1))
		elif r3:
			print(r3.group(1))



