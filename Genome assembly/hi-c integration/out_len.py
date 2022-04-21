#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # fasta

def read_fasta_len(fasta):
	scaf = {}     # to calculate the length of each scaffold
	name = ''
	seq= ''
	a=0
	for line in fasta:
	    line=line.strip()
	    if '>' in line:
	        scaf[name] = len(seq)
	        seq = ''
	        name = line[1:]
	        a = 1
	        continue
	    elif a and '>' not in line:
	        seq += line
	scaf[name] = len(seq)
	scaf.pop('')
	for k in scaf:
	    print('>'+k,'\t',scaf[k])

read_fasta_len(i1)

