#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 
o1 = open(sys.argv[2], 'w')  # go
o2 = open(sys.argv[3], 'w')  # ipr

print('GID\tGO\tEVIDENCE')
for line in i1:
	if '#Chr' not in line:
		line = line.strip().split('\t')
		go = line[4]
		ipr = line[5]
		if go != '-':	
			o1.write('%s\t' % line[3])
			for k in range(0, len(go.strip().split(';'))):
				if go.strip().split(';')[k] != '':
					print(line[3], end='\t')
					print(go.strip().split(';')[k].split()[0], "IEA", sep='\t')
					if k != (len(go.strip().split(';')) - 2):
						o1.write('%s,' % go.strip().split(';')[k].split()[0])
					else:
						o1.write('%s\n' % go.strip().split(';')[k].split()[0])
		if ipr != '-':
			o2.write('%s\t' % line[3])
			for k in range(0, len(ipr.split(';'))):
				if ipr.strip().split(';')[k] != '':
					if k != (len(ipr.split(';')) - 2):
						o2.write('%s\t' % ipr.split(';')[k].split()[0])
					else:
						o2.write('%s\n' % ipr.split(';')[k].split()[0])

