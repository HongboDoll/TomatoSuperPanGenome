#!/public/software/bin/python
# extract each seq

import sys

f1=open(sys.argv[1])
fa=sys.argv[2]
f2=open(sys.argv[3],'w')

number=0
for eachline in f1:
	line=eachline.split()
	name=line[0][1:]
	number+=1
	f2.write("%s%s %s %s%s%s\n"%("/vol1/agis/huangsanwen_group/lihongbo/bin/give_me_one_seq.pl ",fa,name," >./seq/",number,".fa"))


f1.close()
f2.close() 
