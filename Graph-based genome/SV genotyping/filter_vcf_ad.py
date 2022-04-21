#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 

for line in i1:
	if '#' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		dp = int(line[-1].split(':')[1])
		ad = line[-1].split(':')[2].split(',')
		if dp >= 2:
			if int(ad[0]) >= 2 or int(ad[1]) >= 2:
				print('\t'.join(line))

