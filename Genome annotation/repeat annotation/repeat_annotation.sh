#!/bin/bash

ref=S_galapagense_canu_pilon.chr.fasta

## 1.RepeatScout  (v1.0.5) 
build_lmer_table   -sequence   $ref  -freq  ${ref}.freq
RepeatScout -sequence  $ref  -output    fa_file.out  -freq   ${ref}.freq


## 2.LTR-FINDER (v1.07)
ltr_finder  -S  6  -tRNA  tRNAdb   -a  ps_scan   >  ${ref}.LTR.alignment

## 3.MITE-hunter： 
MITE_Hunter_manager.pl  -x pcfg  -i  ${ref}   -n 10

## 4.PILER-DF：(v1.05)
pals -self  ${ref}  -out   ${ref}.seq.hit
piler2   -trs  ${ref}.seq.hit    -out   ${ref}.seq.hit.trs
piler2     -trs2fasta    ${ref}.seq.hit.trs   -seq  ${ref}    -path     fams

## 5.PASTEClassifier.py :
PASTEClassifier.py -i lib -C cfg  > PASTEC.log

## 6.RepeatMasker： (v4.0.5)
Repeatmasker  -nolow -no_is -norna -engine wublast  -parallel 5 -qq   -frag 20000    ${ref}  >${ref}.log
