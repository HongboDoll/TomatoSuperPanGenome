#!/public/software/bin/python
#extrat ok gene from snap result


import sys

f1=open(sys.argv[1])
f2=open(sys.argv[2])
f3=open(sys.argv[3],'w')


dict={}

for eachline in f1:
	line=eachline.split()
	key=line[0]
	dict[key]=key
	
for q in f2:
	i=q.split()
#	if ">" in q:
#		f3.write("%s"%(q))
#	else:
#		if i[3] in dict and :
#			f3.write("%s"%(q))
#			number=1
	if i[2]=="CDS":
		m=i[8].split(";")[1][7:]
		if m in dict:
			f3.write("%s"%(q))


f1.close()
f2.close()
f3.close()
