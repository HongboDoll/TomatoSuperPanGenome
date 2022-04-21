# encoding : utf-8

import sys

input=open(sys.argv[1])#gff file

dict={}
list=[]
for eachline in input:
	line=eachline.split()
	if line[2]=='CDS':
		key=line[8].split(";")[1][7:]
#		print(key)
#		key=line[11]
		if key in dict:
			if line[6]=="+":
				w=line[3]+'\t'+line[4]
				dict[key].append(w)
			if line[6]=="-":
				w=line[4]+'\t'+line[3]
				dict[key].append(w)
		if key not in dict:
			dict[key]=[]
			dict[key].append(line[0])
#			dict[key].append(key)
			dict[key].append(line[6])
			if line[6]=="+":
				w=line[3]+'\t'+line[4]
				dict[key].append(w)
				list.append(key)
			if line[6]=="-":
				w=line[4]+'\t'+line[3]
				dict[key].append(w)
				list.append(key)
for i in list:
	if len(dict[i])==3:
		print('>'+dict[i][0])
		w='Esngl'+'\t'+dict[i][2]+'\t'+i
		print(w)
	if len(dict[i])>3:
		print('>'+dict[i][0])
		if dict[i][1]=="+" or dict[i][1]=="-":
			for j in dict[i][2:]:
				if dict[i].index(j)==2:
					w='Einit'+'\t'+j+'\t'+i
					print(w)
				elif dict[i].index(j)==len(dict[i])-1:
					w='Eterm'+'\t'+j+'\t'+i
					print(w)
				else:
					w='Exon'+'\t'+j+'\t'+i
					print(w)
#		if dict[i][1]=="-":
#			dict[i].reverse()
#			for j in dict[i][0:-2]:
#				if dict[i].index(j)==0:
#					w='Eterm'+'\t'+j+'\t'+i
#					print(w)
#				elif dict[i].index(j)==len(dict[i])-3:
#					w='Einit'+'\t'+j+'\t'+i
#					print(w)
#				else:
#					w='Exon'+'\t'+j+'\t'+i
#					print(w)
input.close()
