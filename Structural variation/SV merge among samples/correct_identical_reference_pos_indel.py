#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 12sol.ZY56.Filter.Coline.indel.merge.vcf
i2 = open(sys.argv[2]) # 

ll = ''
for line in i2:
	ll += line

d = []
for line in i1:
	line = line.strip()
	sline = line.split()
	if ll.strip().split()[0] == sline[0] and int(ll.strip().split()[1]) == int(sline[1]):
		d.append(line)
ref_allel = ''
len_ref = []
for k in d:
	len_ref.append(len(k.split()[3]))

ref_allel = d[len_ref.index(max(len_ref))].split()[3]
gt = ['0/0' for i in range(0, len(d[0].split()[9:]))]

chrom = d[0].split()[0]
pos = d[0].split()[1]

gt_order = 1
gt_order_output = 1
alt = ''
for kk in d:
	k = kk.split()
	if k[3] == ref_allel:
		alt += k[4]+','
		if len(k[4].split(',')) == 1:
			for j in range(0, len(k[9:])):
				if k[j+9] == str(gt_order)+'/'+str(gt_order):
					gt[j] = str(gt_order_output)+'/'+str(gt_order_output)
			gt_order_output += 1
		else:
			for i in range(0, len(k[4].split(','))):
				for j in range(0, len(k[9:])):
					if k[j+9] == str(gt_order+i)+'/'+str(gt_order+i):
						gt[j] = str(gt_order_output+i)+'/'+str(gt_order_output+i)
			gt_order_output += len(k[4].split(','))
	elif k[3] != ref_allel:
		if len(k[4].split(',')) == 1:
			if len(k[4]) <= len(k[3]):
				alt += k[4]+ref_allel[len(k[3]):]+','
				for j in range(0, len(k[9:])):
					if k[j+9] == str(gt_order)+'/'+str(gt_order):
						gt[j] = str(gt_order_output)+'/'+str(gt_order_output)
				gt_order_output += 1
			elif len(k[4]) > len(k[3]):
				alt += k[4]+ref_allel[len(k[3]):]+','
				for j in range(0, len(k[9:])):
					if k[j+9] == str(gt_order)+'/'+str(gt_order):
						gt[j] = str(gt_order_output)+'/'+str(gt_order_output)
				gt_order_output += 1
		else:
			for i in range(0, len(k[4].split(','))):
				if len(k[4].split(',')[i]) <= len(k[3]):
					alt += k[4].split(',')[i]+ref_allel[len(k[3]):]+','
					for j in range(0, len(k[9:])):
						if k[j+9] == str(gt_order+i)+'/'+str(gt_order+i):
							gt[j] = str(gt_order_output+i)+'/'+str(gt_order_output+i)
					
				elif len(k[4].split(',')[i]) > len(k[3]):
					alt += k[4].split(',')[i]+ref_allel[len(k[3]):]+','
					for j in range(0, len(k[9:])):
						if k[j+9] == str(gt_order+i)+'/'+str(gt_order+i):
							gt[j] = str(gt_order_output+i)+'/'+str(gt_order_output+i)

			gt_order_output += len(k[4].split(','))
	
alt = alt[:-1]

print(chrom, pos, '.', ref_allel, alt, '30', '.', 'DP=20', 'GT', sep='\t', end='\t')
for k in gt:
	print(k, end='\t')
print('')

	
