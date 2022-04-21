#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # M82_ZY56syri.vcf

for line in i1:
	if '#' not in line and line.split()[4] == '<TRANS>':
		line = line.strip().split()
		ref_end = line[-1].split(';')[0].split('=')[1]
		alt_chr = line[-1].split(';')[1].split('=')[1]
		alt_start = line[-1].split(';')[2].split('=')[1]
		alt_end = line[-1].split(';')[3].split('=')[1]
		print('#', line[0], line[1], ref_end, '-', alt_chr, alt_start, alt_end, sep='\t')
