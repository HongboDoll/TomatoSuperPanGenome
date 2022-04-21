#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 362_annotated_metabolite_1MbWindow_count.xls
i2 = open(sys.argv[2])  # 362_annotated_metabolite_ID_name_corres.xls_4_plot

corres = {}
for line in i2:
	line = line.strip().split('\t')
	corres[line[0]] = line[2]

print('CHR\tSTART\tEND\tOthers\tAA and NA derivative\tGlycoalkaloid\tPolyamine\tFlavonoid\tPolyphenol\tPhytohormone\tSaccharides')
name = {}
for line in i1:
	if 'CHR' in line:
		line = line.strip().split()
		for k in range(2, len(line)):
			name[k] = line[k]
	else:
		ot = 0
		aa = 0
		gly = 0
		polya = 0
		fla = 0
		polyp = 0
		phy = 0
		sac = 0
		line = line.strip().split()
		print(line[0], line[1].split('~')[0], line[1].split('~')[1], sep='\t', end='\t')
		for k in range(2, len(line)):
			if corres[name[k]] == "Others":
				ot += int(line[k])
			elif corres[name[k]] == "AA and NA derivative":
				aa += int(line[k])
			elif corres[name[k]] == "Glycoalkaloid":
				gly += int(line[k])
			elif corres[name[k]] == "Polyamine":
				polya += int(line[k])
			elif corres[name[k]] == "Flavonoid":
				fla += int(line[k])
			elif corres[name[k]] == "Polyphenol":
				polyp += int(line[k])
			elif corres[name[k]] == "Phytohormone":
				phy += int(line[k])
			elif corres[name[k]] == "Saccharides":
				sac += int(line[k])
		print(ot, aa, gly, polya, fla, polyp, phy, sac, sep='\t')

