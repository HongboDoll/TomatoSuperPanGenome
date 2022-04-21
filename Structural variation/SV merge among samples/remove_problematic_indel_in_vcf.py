#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 12sol.ZY56.Filter.Coline.indel.problem.pos.xls
i2 = open(sys.argv[2])  # 12sol.ZY56.Filter.Coline.indel.merge.vcf

d = {}
for line in i1:
	line = line.strip().split()
	d[line[0]+'~'+line[1]] = ''

for line in i2:
	line = line.strip()
	sline = line.split()
	if sline[0]+'~'+sline[1] not in d:
		print(line)
