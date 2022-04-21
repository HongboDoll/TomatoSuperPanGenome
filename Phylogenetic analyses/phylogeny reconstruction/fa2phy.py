#!/usr/bin/env python

import re
import argparse
parser = argparse.ArgumentParser("Convert fasta format to phylip format")
parser.add_argument("-i", type = str, help = "Input the name of the fasta file")
parser.add_argument("-o", type = str, help = "Output the name of the phylip file")
args = parser.parse_args()

with open(args.i, 'r') as fin:
    sequences = [(m.group(1), ''.join(m.group(2).split()))
    for m in re.finditer(r'(?m)^>([^ \n]+)[^\n]*([^>]*)', fin.read())]
with open(args.o, 'w') as fout:
    fout.write('%d %d\n' % (len(sequences), len(sequences[0][1])))
    for item in sequences:
        fout.write('%-40s %s\n' % item)

