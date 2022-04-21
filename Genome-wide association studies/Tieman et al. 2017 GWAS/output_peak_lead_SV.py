#!/usr/bin/env python3

import sys

i1 = sys.stdin  # awk '$4<2.69E-5' SlFM1804.PAV.out.txt | 
i2 = int(sys.argv[1])  # 1000000 1-Mb window size
o1 = open(sys.argv[2], 'w') # peak signal with no redandancy (1-Mb merged)
o2 = open(sys.argv[3], 'w') # number of signals per 1-Mb window of each chromosome


chr_len = {1:97000000,2:56000000,3:67000000,4:66000000,5:69000000,6:52000000,7:71000000,8:66000000,9:72000000,10:69000000,11:55000000,12:68000000}

chr_window = {}
for k in range(1,13):
	chr_window[k] = {}

for k in chr_len:
	for n in range(0, int(chr_len[k]/i2)):
		start = (n*i2)+1
		end = (n+1)*i2
		chr_window[k][str(start)+'~'+str(end)] = {}
	
for line in i1:
	line = line.strip().split()
	for k in chr_window[int(line[1])]:
		if int(k.split('~')[0]) <= int(line[2]) and int(k.split('~')[1]) >= int(line[2]):
			chr_window[int(line[1])][k][line[0]+'~'+line[1]+'~'+line[2]+'~'+line[3]] = ''

for k in chr_window:
	for j in chr_window[k]:
		l = sorted(list(chr_window[k][j].keys()), key = lambda i: (float(i.split('~')[3])), reverse = False)
		if l and len(l) >=1:
			o1.write('%s\n' % '\t'.join(l[0].split('~')))
			o2.write('%s\t%s\t1\n' % (k, j))
		else:
			o2.write('%s\t%s\t0\n' % (k, j))

