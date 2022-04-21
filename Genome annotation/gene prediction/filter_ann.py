#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # genome.zff2keep
i2 = open(sys.argv[2])  # genome.ann

d = {}
for line in i1:
	d[line.strip().split()[1]] = ''

for line in i2:
	if '>' in line:
		print(line.strip())
	else:
		line = line.strip()
		sline = line.split()
		if sline[-1] in d:
			print(line)
