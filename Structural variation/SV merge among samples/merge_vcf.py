#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # existed vcf

d1 = {}
spe = ['Heinz','M82','ZY57','ZY58','ZY60','ZY61','ZY62','ZY65','ZY59','Spenn','ZY63','ZY64']
for line in i1:
	line = line.strip().split()
	if line[0]+'~'+line[1]+'~'+line[3] not in d1:
		d1[line[0]+'~'+line[1]+'~'+line[3]] = []
		d1[line[0]+'~'+line[1]+'~'+line[3]].append(line[4]+'~'+line[-1])
	else:
		d1[line[0]+'~'+line[1]+'~'+line[3]].append(line[4]+'~'+line[-1])

print('1CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHeinz\tM82\tZY57\tZY58\tZY60\tZY61\tZY62\tZY65\tZY59\tSpenn\tZY63\tZY64')
for k in d1:
	l = []
	geno = []
	for j in d1[k]:
		l.append(j.split('~')[1])
		geno.append(j.split('~')[0])
	geno = sorted(list(set(geno)))
	genoseq =''
	for kk in geno:
		genoseq += kk+','
	genoseq = genoseq[:-1]
	ll = set(spe) - set(l)
	for i in ll:
		d1[k].append(k.split('~')[2]+'~'+i)
	print(k.split('~')[0],k.split('~')[1],'.',k.split('~')[2],genoseq,'30','.','DP=100','GT',sep='\t',end='\t')
	for s in spe:
		for i in d1[k]:
			if s == i.split('~')[1]:
				if i.split('~')[0] == k.split('~')[2]:
					print('0/0', end='\t')
				else:
					num = str(geno.index(i.split('~')[0]) + 1)
					print(num+'/'+num, end='\t')
	print('')

					
					
					
