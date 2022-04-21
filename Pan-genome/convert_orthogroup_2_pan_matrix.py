#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # Orthogroups.tsv

sample = {}
#nn = 1
for line in i1:
	if 'Orthogroup' in line:
		line = line.strip().split('\t')
		print('\t'.join(line).replace('Orthogroup', 'Cluster'))
		for n in range(1, len(line)):
			sample[n] = line[n]
	else:
		line = line.split('\t')
		print(line[0],end='\t')
#		print(nn, end='\t')
#		nn += 1
		#if line[45] == '\n':
	#		print('11,',line[45].strip())
		for n in range(1, len(line)):
	#		pass
			#if line[n] == '':
			#	print(n, line[n])
			if line[n] == '':
				if n != len(line)-1:
					print('-',end='\t')
				else:
					print('-')
			elif line[n] != '\n' and line[n] != '':
				if n != len(line)-1:
					print(line[n].strip().replace(' ', ''),end='\t')
				else:
					print(line[n].strip().replace(' ', ''))
			elif line[n] == '\n':
				print('-')
				
