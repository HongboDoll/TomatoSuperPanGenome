#!/bin/bash

ls emmax_snpgwas/*txt | while read i
do
	name=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
	awk '$4<7.58E-7' $i | ./output_peak_lead_SV.py 1000000 emmax_snpgwas/${name}.GAL.1MbNonredundant.xls emmax_snpgwas/${name}.GAL.1MbWindow.count
done
#
ls sv_results/*txt | while read i
do
	name=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
	awk '$4<3.19E-5' $i | ./output_peak_lead_SV.py 1000000 sv_results/${name}.PAV.1MbNonredundant.xls sv_results/${name}.PAV.1MbWindow.count
done

ls emmax_snpgwas/*txt | while read i
do
	name=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
	./find_sv_snpGWAS_overlap_uniq_signals_emmax.py sv_results/${name}*1MbNonredundant.xls emmax_snpgwas/${name}*1MbNonredundant.xls > sv_results/${name}_share_uniq.xls
done
