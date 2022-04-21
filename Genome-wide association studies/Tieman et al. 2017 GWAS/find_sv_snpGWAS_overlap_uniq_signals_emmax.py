#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # svgwas_signals.xls
i2 = open(sys.argv[2])  # snpgwas_signals.xls

sv = {}
for line in i1:
	line = line.strip().split()
	sv[line[1].replace('\"','')+'_'+line[2]] = line

snp = {}
for line in i2:
	line = line.strip().split()
	snp[line[1].replace('\"','')+'_'+line[2]] = line

for k in sv:
	chr_name = k.split('_')[0]
	pos = int(k.split('_')[1])
	shared_flag = 0
	for j in snp:
		chr_name2 = j.split('_')[0]
		pos2 = int(j.split('_')[1])
		if chr_name == chr_name2:
			if pos - 500000 <= pos2 and pos + 500000 >= pos2:
				shared_flag += 1
			elif pos - 1000000 <= pos2 and pos >= pos2:
				shared_flag += 1
			elif pos <= pos2 and pos + 1000000 >= pos2:
				shared_flag += 1
	if shared_flag:
		print(','.join(sv[k]), 'shared', sep='\t')
	else:
		print(','.join(sv[k]), 'svgwas_unique', sep='\t')

