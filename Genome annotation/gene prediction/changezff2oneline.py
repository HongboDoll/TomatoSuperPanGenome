#change the format of zff to one name



import sys

f1=open(sys.argv[1])
f2=open(sys.argv[2])
f3=open(sys.argv[3],'w')

dict={}
for eachline in f1:
	line=eachline.split()
	key=line[0]
	dict[key]=key

for p in f2:
	i=p.split()
	#print(len(p))
	if ">" in p and i[0] in dict:
		f3.write("%s"%(p))
		del(dict[i[0]])
	elif ">" not in p:
		f3.write("%s"%(p))

f1.close()
f2.close()
f3.close()
