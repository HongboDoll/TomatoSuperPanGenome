#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 
o1 = open(sys.argv[2], 'w') # stat

length = 0
length_95 = 0
flag = 0
for line in i1:
	if '#' in line:
		flag += 1
	else:
		flag = 0
		line = line.strip().split()
#		if abs(int(line[0]) - int(line[-1].split(':')[1].split('-')[0])) <= 10000000 and abs(int(line[1]) - int(line[-1].split(':')[1].split('-')[1])) <= 10000000:
		if float(line[6]) >= 95:
			length += int(line[5])
		else:
			length_95 += int(line[5])
	
if flag:
	o1.write('Unaligned\n')
else:
	o1.write('%s %s\n' % (length, length_95))

if length >= 301:
	print('+')
else:
	print('-')



