#!/public/software/bin/python

import sys
import copy
import re

fin=open(sys.argv[1])
fout=open(sys.argv[2],'w')

for i in fin:
	if '#' not in i:
		a=i.split('\t')
		if a[2]=='mRNA':
			b=copy.deepcopy(a)
			a[2]='gene'
#			print a
			fout.write('\t'.join(a))
			blist=b[-1].split(';')
			mrna_left=re.sub('gene','m',blist[0])
			mrna_right=re.sub('Name','Parent',blist[1])
			tmp_str=mrna_left+';'+mrna_right
			b=b[:8]+[tmp_str]
#			print b
			fout.write('\t'.join(b))
		else:
			c=copy.deepcopy(a)
			c[2]='exon'
			clist=c[-1].split(';')
			left=re.sub('cds','exon',clist[0])
			right='Parent='+mrna_left[3:]
			tmp_str=left+';'+right
			c=c[:8]+[tmp_str]
	#		print c
			fout.write('\t'.join(c)+'\n')
			a=a[:8]+[a[-1].split(';')[0]+';'+right]
			fout.write('\t'.join(a)+'\n')
#			print a
	
