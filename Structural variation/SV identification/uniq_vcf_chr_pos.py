#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # existed vcf

d = {}
for line in i1:
	line = line.strip()
	sline = line.split()
	if sline[0]+'_'+sline[1] not in d:
		print(line)
	d[sline[0]+'_'+sline[1]] = ''
										
											
