#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 

def read_fasta(fasta):
	scaf = {}     # to calculate the length of each scaffold
	name = ''
	seq= ''
	a=0
	for line in fasta:
	    line=line.strip()
	    if '>' in line:
	        scaf[name] = seq
	        seq = ''
	        name = line[1:]
	        a = 1
	        continue
	    elif a and '>' not in line:
	        seq += line
	scaf[name] = seq
	scaf.pop('')
	return scaf

scaf = read_fasta(i1)

for k in scaf:
	if 'TRINITY_GG' not in k:
		print(">"+k)
		print(scaf[k])
