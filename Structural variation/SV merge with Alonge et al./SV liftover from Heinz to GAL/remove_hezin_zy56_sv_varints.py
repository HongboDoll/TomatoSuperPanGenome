#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # heinz_zy56_sv.xls
i2 = open(sys.argv[2])  # merged.ont.v1.0.convert_GAL_filter300.vcf
o = open(sys.argv[3], 'w') # output stat

d = {}
for line in i1:
	line = line.strip().split()
	if line[0] not in d:
		d[line[0]] = {}
		d[line[0]][int(line[1])] = ''
	else:
		d[line[0]][int(line[1])] = ''

		
for line in i2:
	if '#' in line:
		print(line.strip())
	else:
		line = line.strip().split()
		f = 0
		if line[0] in d:
			for k in d[line[0]]:
				if int(line[1]) - 20 <= k and int(line[1]) + 20 >= k and not f:
					f += 1
					o.write('%s\t%s\t%s\n' % (line[0], line[1], line[2]))
		if not f:
			print('\t'.join(line))

