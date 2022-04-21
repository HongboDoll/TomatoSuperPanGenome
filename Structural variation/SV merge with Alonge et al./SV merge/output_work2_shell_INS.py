#!/usr/bin/env python3

import sys

n = int(sys.argv[1])  # 
n1 = (n-1)*580+1
n2 = n*580

print("""#!/bin/bash

grep -v '#' 12_100_merged_ZY56_INS.vcf | """, end='')

print("sed -n '"+str(n1)+','+str(n2)+"p' |", end='')

print(""" while read i
do
		chr=`echo $i | awk '{print $1}'`
		pos=`echo $i | awk '{print $2}'`
		id=`echo $i | awk '{print $3}'`
		ref=`echo $i | awk '{print $4}'`
		length=`echo $ref | wc | awk '{print $3}'`
		len=`expr $length - 2`
		end=`expr $pos + $len`
		ref_seq=`samtools faidx gal_ref/${chr}.fa ${chr}:${pos}-${end} | grep -v '>'  | sed ':a;N;s/\\n//g;ta'`
		paste <(echo $i | cut -d ' ' -f '1-3') <(echo $ref_seq) <(echo $i | cut -d ' ' -f '5-999') | sed 's/\ /\\t/g' > revise_ref/${chr}_${pos}_${id}_INS
done
""")
