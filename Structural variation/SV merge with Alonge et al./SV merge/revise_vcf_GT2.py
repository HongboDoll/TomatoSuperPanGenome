#!/usr/bin/env python3
import sys

i1 = open(sys.argv[1])  # merged.vcf without genotype
i2 = open(sys.argv[2])  # 12sol original.vcf1 with genotype 
i3 = open(sys.argv[3])  # 100cell original.vcf2 with genotype
i4 = open(sys.argv[4])  # head_12
i5 = open(sys.argv[5])  # head_100

acc1 = []
acc2 = []
for line in i4:
	acc1 = line.strip().split()

for line in i5:
	acc2 = line.strip().split()

geno1 = {}
for line in i2:
	if '#' not in line:
		line = line.strip().split()
		if line[0]+'_'+line[1] not in geno1:
			geno1[line[0]+'_'+line[1]] = line[9:]
		else:
			for n in range(0, len(geno1[line[0]+'_'+line[1]])):
				if line[n+9] != '0/0' and line[n+9] != './.' and (geno1[line[0]+'_'+line[1]][n] == '0/0' or geno1[line[0]+'_'+line[1]][n] == './.'):
					geno1[line[0]+'_'+line[1]][n] = line[n+9]
		for k in geno1[line[0]+'_'+line[1]]:
			if k != '0/0' and k != './.' and k != '0/1' and k != '1/0':
				k = '1/1'

geno2 = {}
for line in i3:
	if '#' not in line:
		line = line.strip().split()
		if line[0]+'_'+line[1] not in geno2:
			geno2[line[0]+'_'+line[1]] = line[9:]
		else:
			for n in range(0, len(geno2[line[0]+'_'+line[1]])):
				if line[n+9] != '0/0' and line[n+9] != './.' and (geno2[line[0]+'_'+line[1]][n] == '0/0' or geno2[line[0]+'_'+line[1]][n] == './.'):
					geno2[line[0]+'_'+line[1]][n] = line[n+9]
		for k in geno2[line[0]+'_'+line[1]]:
			if k != '0/0' and k != './.' and k != '0/1' and k != '1/0':
				k = '1/1'
			
for line in i1:
	if '#' in line and '#CHROM' not in line:
		print(line.strip())
	elif '#' in line and '#CHROM' in line:
		line = line.strip().split()
		print('\t'.join(line[0:9]), end='\t')
		for k in acc1:
			print(k, end='\t')
		for n in range(0, len(acc2)):
			if n != (len(acc2) - 1):
				print(acc2[n], end='\t')
			else:
				print(acc2[n])
	else:
		line = line.strip().split()
		gt = []
		if 'SUPP=1' in line[-4]:
			if line[-2].split(':')[-1] == 'NAN':
				for n in range(0, 12):
					gt.append('0/0')
				for n in geno2[line[-1].split(':')[-1].split('-')[0]]:
					gt.append(n)
			else:
				for n in geno1[line[-2].split(':')[-1].split('-')[0]]:
					gt.append(n)
				for n in range(0, 100):
					gt.append('0/0')
		elif 'SUPP=2' in line[-4]:
			for n in geno1[line[-2].split(':')[-1].split('-')[0]]:
				gt.append(n)
			for n in geno2[line[-1].split(':')[-1].split('-')[0]]:
				gt.append(n)

		print('\t'.join(line[0:6]), '.', '.', 'GT', '\t'.join(gt), sep='\t')

