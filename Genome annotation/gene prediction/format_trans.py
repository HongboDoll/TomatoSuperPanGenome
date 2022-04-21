#!/usr/bin/python

import sys

fin=open(sys.argv[1])
fout=open(sys.argv[2],'w')

for i in fin:
	if '>' in i:
		Chr=i.strip()[1:]
	else:
		
		if 'Esngl' not in i and 'Eterm' not in i:
			a=i.split()
			item='%s\t%s\t%s\n'%(Chr,a[1],a[2])
			fout.write(item)
		else:
			a=i.split()
			item='%s\t%s\t%s\n\n'%(Chr,a[1],a[2])
			fout.write(item)
			
		

