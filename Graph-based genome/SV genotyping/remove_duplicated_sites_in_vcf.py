#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # duplicated pos
i2 = open(sys.argv[2])  # vcf

pos = {}
pos_count = {}
for line in i1:
	line = line.strip().split()
	pos[line[1]+'~'+line[2]] = 0
	pos_count[line[1]+'~'+line[2]] = int(line[0])

l = []
for line in i2:
	if '#' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		if line[0]+'~'+line[1] in pos:
			if pos[line[0]+'~'+line[1]] == 0:
				pos[line[0]+'~'+line[1]] += 1
				print('\t'.join(line[0:9]), end='\t')
				for n in range(9, len(line)):
					l.append(line[n])
			elif pos[line[0]+'~'+line[1]] != 0 and pos[line[0]+'~'+line[1]] != (pos_count[line[0]+'~'+line[1]] - 1):
				pos[line[0]+'~'+line[1]] += 1
				for n in range(9, len(line)):
					if line[n].split(':')[0] != './.':
#						print('')
#						print(n, line[n], l[(n-9)].split(':')[0])
						if '1' not in l[(n-9)].split(':')[0]:
							l[(n-9)] = line[n]
#				print(l)
			elif pos[line[0]+'~'+line[1]] == (pos_count[line[0]+'~'+line[1]] - 1):
				for n in range(9, len(line)):
					if line[n].split(':')[0] != './.':
#						print('')
#						print(n, line[n])
						if '1' not in l[(n-9)].split(':')[0]:
							l[(n-9)] = line[n]
#				print(l)
				for k in l:
					print(k, end='\t')
				print('')
				l = []
		else:
			print('\t'.join(line))

