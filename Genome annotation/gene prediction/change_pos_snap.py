#!/usr/bin/python
#change gene pos to fit the zff pos 
#########2016.11.13#######################


import sys

f1=open(sys.argv[1]) # gene.fa 
f2=open(sys.argv[2]) # pos_uniq
f3=open(sys.argv[3],'w')

scaf = {}     
name = ''
seq= ''
for line in f1:
    line=line.strip()
    if '>' in line:
        scaf[name] = seq
        seq = ''
        name = line.split()[0]
        a = 1
        continue
    elif a and '>' not in line:
        seq += line
scaf[name] = seq
scaf.pop('')

for eachline in f2:
	i=eachline.split()
	if ">" in eachline:
#		print(i[0])
		if i[0] in scaf:
			f3.write('%s\n%s\n'%(i[0],scaf[i[0]]))	
		

f1.close()
f2.close()
		
		
		
	
