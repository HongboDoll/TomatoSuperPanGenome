#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # cm.ZY65_ZY56_svmu.txt
i2 = open(sys.argv[2])  # ZY65_ZY56.delta

syn_block = {}
syn_uniq = {}
for line in i1:
	line = line.strip().split()
	syn_block[line[0]+'~'+line[1]+'~'+line[2]+'~'+line[3]+'~'+line[4]+'~'+line[5]] = ''
	syn_uniq[line[0]+'~'+line[3]] = ''

d = {}
name = ''
d[name] = {}
seq = ''
pos = []
a = 0
corre1 = {}
corre2 = {}
for line in i2:
	if '/' in line or 'NUCMER' in line:
		print(line.strip())
	else:
		line = line.strip()
		if not a and '>' in line:
			corre1[line.split()[0][1:]+'~'+line.split()[1]] = line
			d[name][seq] = pos
			d[line] = {}
			seq = ''
			pos = []
			name = line
			a = 1
		elif '>' not in line and len(line.split()) == 7:	
			corre2['~'.join(line.split()[0:4])] = line
			d[name][seq] = pos
			seq = line
			d[name][seq] = []
			pos = []
		elif '>' not in line and len(line.split()) == 1:
			pos.append(line)
		elif a and '>' in line:
			corre1[line.split()[0][1:]+'~'+line.split()[1]] = line
			d[name][seq] = pos
			d[line] = {}
			seq = ''
			pos = []
			name = line
d[name][seq] = pos
d.pop('')
for i in syn_uniq:
	if corre1[i.split('~')[0]+'~'+i.split('~')[1]] in d:
		print(corre1[i.split('~')[0]+'~'+i.split('~')[1]])
		for k in syn_block:
			if k.split('~')[0] == i.split('~')[0] and k.split('~')[3] == i.split('~')[1]:
				if corre2[k.split('~')[1]+'~'+k.split('~')[2]+'~'+k.split('~')[4]+'~'+k.split('~')[5]] in d[corre1[k.split('~')[0]+'~'+k.split('~')[3]]]:
					print(corre2[k.split('~')[1]+'~'+k.split('~')[2]+'~'+k.split('~')[4]+'~'+k.split('~')[5]])
				for j in d[corre1[k.split('~')[0]+'~'+k.split('~')[3]]][corre2[k.split('~')[1]+'~'+k.split('~')[2]+'~'+k.split('~')[4]+'~'+k.split('~')[5]]]:
					print(j)
		
