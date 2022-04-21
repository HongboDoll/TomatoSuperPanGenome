#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code take output from show-snps.
The options should be set as:
show-snps -Clr -x 1  -T mum.delta.filter  >mum.delta.filterX.snps
Usage:
python3.4 MUMmerSNPs2VCF.py mum.delta.filterX.snps mum_filterX.snps.vcf


Keywords: MUMmer show-snps VCF MUMmer2VCF
"""


from sys import argv
script, input_file, output_file = argv

OUT = open(output_file,"w")
OUT.write('##fileformat=VCFv4.1'+"\n")
OUT.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+"\n")
out_line = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample1']
OUT.write("\t".join(out_line)+"\n")						
vcf_out = []


def check_buff(indel_buff_in):
	allele_ref = indel_buff_in[0][1]
	allele_alt = indel_buff_in[0][2]
	ref_id  = indel_buff_in[0][12]	
	if allele_ref == '.':
		# insertion 
		pos  = indel_buff_in[0][0]
		# In MUMmer format, the coordinate of '.' is the coordinate of the last nt so, this position is kept.
		ref_start = indel_buff_in[0][8][0]
		direction = indel_buff_in[0][11]
		alt_out = ''
		if direction == '-1':
			for line_l in indel_buff_in :
				alt_out = line_l[2]+alt_out
		else :
			for line_l in indel_buff_in :
				alt_out += line_l[2]
		alt_out = ref_start+alt_out
		out_line = [ref_id,pos,'.',ref_start,alt_out,'30','PASS','.','GT','1/1']
		vcf_out.append(out_line)
	elif allele_alt == '.':
		# deletion
		pos  = str(int(indel_buff_in[0][0])-1)
		# the coordinate here in the reference is correct, but we need the coordinate of last nt.
		# In VCF format, we need check the last nt.
		alt_start = indel_buff_in[0][8][0] # first nt in context
		ref_out = alt_start
		for line_l in indel_buff_in :
			ref_out += line_l[1]
		out_line = [ref_id,pos,'.',ref_out,alt_start,'30','PASS','.','GT','1/1']
		vcf_out.append(out_line)
	else :
		sys.exit("Both in and del\n")

#####################################
# initiation 
start  = 0
last_pos = 0
last_ref = ''
in_del_start = 0
indel_buff = []
##################################

with open (input_file,"r") as INPUT:
	for line in INPUT:
		line = line.rstrip()
		if len(line)< 1:
			continue
		elif start == 0 and line[0] == '[':
			start = 1
		elif start == 1:
			line_list = line.split("\t")
			ref_id  = line_list[12]
			pos  = line_list[0]
			allele_ref = line_list[1]
			allele_alt = line_list[2]
			if allele_ref == '.' or allele_alt == '.':
				# insertion	 deletion
				if in_del_start == 0:
					in_del_start = 1
					indel_buff.append(line_list)
				else :
					if allele_ref == '.' :
						if ref_id == last_ref and int(pos) == last_pos :
							indel_buff.append(line_list)
						else : # new insertion
							check_buff(indel_buff)
							indel_buff = []
							indel_buff.append(line_list)
					elif allele_alt == '.':
						if ref_id == last_ref and int(pos) == last_pos + 1:
							indel_buff.append(line_list)
						else:  # new deletion
							check_buff(indel_buff)
							indel_buff = []
							indel_buff.append(line_list)
			else :
				# SNP
				if in_del_start == 1:
					check_buff(indel_buff)
					indel_buff = []
					in_del_start = 0
				## write SNP regard less of last records
				out_line = [ref_id,pos,'.',allele_ref,allele_alt,'30','PASS','.','GT','1/1']
				vcf_out.append(out_line)
			##
			last_pos = int(pos)
			last_ref = ref_id
###############


#  Write VCF
new_list1 = sorted(vcf_out, key=lambda x: int(x[1]))
new_list = sorted(new_list1, key=lambda x: x[0])
for line_new in new_list:
	OUT.write("\t".join(line_new)+"\n")
OUT.close()
## Author : lxue@uga.edu
   
